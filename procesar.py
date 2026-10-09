#!/usr/bin/env python3
"""
procesar.py — Genera el mapa de nubes regional que consume la web app "Cielo".

Se ejecuta en GitHub Actions cada 15 min. Para la región configurada (por defecto
Chile + Patagonia) produce una rejilla regular lat/lon de 0.02° (~2 km) donde cada
celda lleva un código:

    0  sin dato
    1  despejado
    2  despejado dudoso: probabilidad de nube 0.35–0.5 (bruma / cirro tenue posible)
    3  nube baja   (tope < 2.5 km)
    4  nube media  (2.5–6 km, o altura desconocida)
    5  nube alta   (> 6 km)
    6  niebla / estrato nocturno (test BT11.2 − BT3.9; la máscara lo daba claro)

Las nubes se dibujan en su POSICIÓN REAL: cada píxel nuboso se desplaza hacia el
punto subsatélite h·tan(VZA) usando la altura de tope ACHA (corrección de paralaje).

Fuente: GOES-19 ABI (bucket público noaa-goes19 en AWS, sin claves):
    ABI-L2-ACMF   máscara de nubes (BCM binaria + ACM 4 niveles)
    ABI-L2-CMIPF  C14 (11.2 µm) y C07 (3.9 µm)
    ABI-L2-ACHAF  altura del tope nuboso

Salida (directorio --salida):
    data/latest.json            metadatos + lista de cuadros (histórico ~4 h)
    data/f/AAAAMMDDhhmm.bin.gz  rejilla uint8 (filas N→S, columnas O→E), gzip

El histórico se recupera del propio sitio publicado (--sitio), así no se guarda
nada en git y el repositorio no crece.
"""
from __future__ import annotations

import argparse
import gzip
import itertools
import json
import math
import os
import re
import sys
from datetime import datetime, timedelta, timezone

import numpy as np
import xarray as xr
from pyproj import Geod, Proj

GEOD = Geod(ellps="WGS84")
R_TIERRA_KM = 6371.0

REGION = dict(lat_max=-17.5, lat_min=-56.0, lon_min=-76.0, lon_max=-64.0, paso=0.02)
HISTORICO_HORAS = 4.0
MAX_CUADROS = 20

UMBRAL_NIEBLA_BTD = 2.5     # K
BT11_MIN_NIEBLA = 253.0     # K, evita superficies muy frías (nieve, hielo)
H_BAJA = 2500.0             # m
H_ALTA = 6000.0             # m
VZA_MAX = 72.0              # grados; más allá se marca sin dato
PROB_DUDOSA = 0.35          # prob. de nube a partir de la cual un píxel "despejado" se marca dudoso

CODIGOS = {0: "sin dato", 1: "despejado", 2: "despejado dudoso", 3: "nube baja",
           4: "nube media", 5: "nube alta", 6: "niebla/estrato nocturno"}

PATRON_T = re.compile(r"_s(\d{4})(\d{3})(\d{2})(\d{2})(\d{2})")


def log(*a):
    print(*a, file=sys.stderr, flush=True)


# ─────────────────────────── geometría ───────────────────────────

def ecef(lat_deg, lon_deg, r_km):
    lat, lon = np.radians(lat_deg), np.radians(lon_deg)
    return np.stack([r_km * np.cos(lat) * np.cos(lon),
                     r_km * np.cos(lat) * np.sin(lon),
                     r_km * np.sin(lat)], axis=-1)


def vza(lat, lon, lon_sat, alt_km=35786.0):
    lat = np.asarray(lat, float)
    lon = np.asarray(lon, float)
    p = ecef(lat, lon, R_TIERRA_KM)
    s = ecef(np.zeros_like(lat), np.full_like(lon, lon_sat), R_TIERRA_KM + alt_km)
    v = s - p
    c = np.sum(v * p, -1) / (np.linalg.norm(v, axis=-1) * np.linalg.norm(p, axis=-1))
    return np.degrees(np.arccos(np.clip(c, -1, 1)))


def desplazamiento_paralaje(lat, lon, h_m, lon_sat):
    """Devuelve (dlat, dlon) a sumar a la posición aparente para obtener la real."""
    if lat.size == 0:
        return np.zeros(0), np.zeros(0)
    z = vza(lat, lon, lon_sat)
    d = (h_m / 1000.0) * np.tan(np.radians(np.minimum(z, 85.0))) * 1000.0
    az, _, _ = GEOD.inv(lon, lat, np.full_like(lon, lon_sat), np.zeros_like(lat))
    lon2, lat2, _ = GEOD.fwd(lon, lat, az, d)
    dlon = (lon2 - lon + 540.0) % 360.0 - 180.0
    return lat2 - lat, dlon


def elevacion_solar(lat, lon, t: datetime):
    lat = np.asarray(lat, float)
    lon = np.asarray(lon, float)
    doy = t.timetuple().tm_yday
    hora = t.hour + t.minute / 60 + t.second / 3600
    g = 2 * math.pi / 365 * (doy - 1 + (hora - 12) / 24)
    decl = (0.006918 - 0.399912 * math.cos(g) + 0.070257 * math.sin(g)
            - 0.006758 * math.cos(2 * g) + 0.000907 * math.sin(2 * g)
            - 0.002697 * math.cos(3 * g) + 0.00148 * math.sin(3 * g))
    eqt = 229.18 * (0.000075 + 0.001868 * math.cos(g) - 0.032077 * math.sin(g)
                    - 0.014615 * math.cos(2 * g) - 0.040849 * math.sin(2 * g))
    tst = hora * 60 + eqt + 4 * lon
    ha = np.radians(tst / 4 - 180)
    la = np.radians(lat)
    cz = np.sin(la) * math.sin(decl) + np.cos(la) * math.cos(decl) * np.cos(ha)
    return 90 - np.degrees(np.arccos(np.clip(cz, -1, 1)))


# ─────────────────────────── acceso a GOES ───────────────────────────

class GOES:
    def __init__(self, raiz="noaa-goes19", local=False):
        if local:
            import fsspec
            self.fs = fsspec.filesystem("file")
        else:
            import s3fs
            self.fs = s3fs.S3FileSystem(anon=True)
        self.raiz = raiz.rstrip("/")

    def listar(self, producto, desde, hasta, filtro=""):
        out, t = [], desde.replace(minute=0, second=0, microsecond=0)
        while t <= hasta:
            ruta = f"{self.raiz}/{producto}/{t:%Y}/{t.timetuple().tm_yday:03d}/{t:%H}"
            try:
                for f in self.fs.ls(ruta):
                    n = f.rsplit("/", 1)[-1]
                    if filtro and filtro not in n:
                        continue
                    m = PATRON_T.search(n)
                    if m:
                        y, d, hh, mm, ss = map(int, m.groups())
                        ts = datetime(y, 1, 1, hh, mm, ss, tzinfo=timezone.utc) + timedelta(days=d - 1)
                        if desde <= ts <= hasta:
                            out.append((ts, f))
            except (FileNotFoundError, OSError):
                pass
            t += timedelta(hours=1)
        return sorted(out)

    def abrir(self, ruta):
        return xr.open_dataset(self.fs.open(ruta, "rb"), engine="h5netcdf")


def mas_cercano(lista, t, tol_min):
    if not lista:
        return None
    ts, f = min(lista, key=lambda x: abs((x[0] - t).total_seconds()))
    return f if abs((ts - t).total_seconds()) <= tol_min * 60 else None


def proyeccion(ds):
    p = ds["goes_imager_projection"].attrs
    H = float(p["perspective_point_height"])
    lon0 = float(p["longitude_of_projection_origin"])
    pj = Proj(proj="geos", h=H, lon_0=lon0, a=float(p["semi_major_axis"]),
              b=float(p["semi_minor_axis"]), sweep=str(p.get("sweep_angle_axis", "x")))
    return pj, H, lon0


def caja_indices(ds, region, margen_px):
    """Índices (slice_y, slice_x) del fixed grid que cubren la región."""
    pj, H, _ = proyeccion(ds)
    n = 200
    lats = np.concatenate([np.linspace(region["lat_min"], region["lat_max"], n)] * 2 +
                          [np.full(n, region["lat_min"]), np.full(n, region["lat_max"])])
    lons = np.concatenate([np.full(n, region["lon_min"]), np.full(n, region["lon_max"])] +
                          [np.linspace(region["lon_min"], region["lon_max"], n)] * 2)
    x, y = pj(lons, lats)
    ok = np.isfinite(x) & np.isfinite(y) & (np.abs(x) < 1e20)
    xs, ys = ds["x"].values, ds["y"].values
    ix = np.array([np.argmin(np.abs(xs - v / H)) for v in x[ok]])
    iy = np.array([np.argmin(np.abs(ys - v / H)) for v in y[ok]])
    sx = slice(max(ix.min() - margen_px, 0), min(ix.max() + margen_px + 1, xs.size))
    sy = slice(max(iy.min() - margen_px, 0), min(iy.max() + margen_px + 1, ys.size))
    return sy, sx


def latlon_caja(ds, sy, sx):
    pj, H, lon0 = proyeccion(ds)
    X, Y = np.meshgrid(ds["x"].values[sx] * H, ds["y"].values[sy] * H)
    lon, lat = pj(X, Y, inverse=True)
    malo = (np.abs(lon) > 360) | (np.abs(lat) > 90) | ~np.isfinite(lon)
    lon = np.where(malo, np.nan, lon)
    lat = np.where(malo, np.nan, lat)
    return lat, lon, lon0


def leer_caja(ds, var, sy, sx):
    return ds[var].isel(y=sy, x=sx).load().values.astype(np.float32)


def remuestrear(ds, var, lat, lon):
    """Lleva una variable de otra rejilla (ACHA, 10 km) a los píxeles dados (vecino más cercano)."""
    pj, H, _ = proyeccion(ds)
    xm, ym = pj(np.nan_to_num(lon), np.nan_to_num(lat))
    ok = np.isfinite(lat) & np.isfinite(xm) & (np.abs(xm) < 1e20)
    xs, ys = ds["x"].values, ds["y"].values
    def idx(eje, v):
        if eje[0] < eje[-1]:
            i = np.searchsorted(eje, v)
        else:
            i = eje.size - np.searchsorted(eje[::-1], v)
        i = np.clip(i, 1, eje.size - 1)
        izq = eje[i - 1]
        der = eje[i]
        return np.where(np.abs(v - izq) < np.abs(v - der), i - 1, i)
    xi = idx(xs, xm / H)
    yi = idx(ys, ym / H)
    out = np.full(lat.shape, np.nan, np.float32)
    if not ok.any():
        return out
    x0, x1 = int(xi[ok].min()), int(xi[ok].max()) + 1
    y0, y1 = int(yi[ok].min()), int(yi[ok].max()) + 1
    bloque = ds[var].isel(x=slice(x0, x1), y=slice(y0, y1)).load().values.astype(np.float32)
    out[ok] = bloque[yi[ok] - y0, xi[ok] - x0]
    return out


# ─────────────────────────── núcleo ───────────────────────────

def procesar_cuadro(goes: GOES, t: datetime, f_acm: str, region=REGION):
    filas = int(round((region["lat_max"] - region["lat_min"]) / region["paso"]))
    cols = int(round((region["lon_max"] - region["lon_min"]) / region["paso"]))

    with goes.abrir(f_acm) as ds:
        sy, sx = caja_indices(ds, region, margen_px=12)
        lat, lon, lon_sat = latlon_caja(ds, sy, sx)
        bcm = leer_caja(ds, "BCM", sy, sx)
        # probabilidad de nube (0–1) del clasificador bayesiano de NOAA
        prob = leer_caja(ds, "Cloud_Probabilities", sy, sx) if "Cloud_Probabilities" in ds \
            else np.full_like(bcm, np.nan)
        x_ref = ds["x"].values[sx]
    log(f"  caja {bcm.shape[0]}×{bcm.shape[1]} píxeles, lon_sat {lon_sat}")

    def canal(nombre):
        lista = goes.listar("ABI-L2-CMIPF", t - timedelta(minutes=12), t + timedelta(minutes=12), f"{nombre}_")
        f = mas_cercano(lista, t, 12)
        if not f:
            log(f"  (sin {nombre})")
            return None
        with goes.abrir(f) as ds:
            if not np.allclose(ds["x"].values[sx], x_ref):
                log(f"  ({nombre}: rejilla distinta, se ignora)")
                return None
            return leer_caja(ds, "CMI", sy, sx)

    bt11 = canal("C14")
    bt39 = canal("C07")

    altura = np.full(lat.shape, np.nan, np.float32)
    f_ht = mas_cercano(goes.listar("ABI-L2-ACHAF", t - timedelta(minutes=70), t + timedelta(minutes=15)), t, 70)
    if f_ht:
        with goes.abrir(f_ht) as ds:
            altura = remuestrear(ds, "HT", lat, lon)
    else:
        log("  (sin ACHA: altura por defecto 4 km)")

    # ── clasificación por píxel ──
    z = vza(np.nan_to_num(lat), np.nan_to_num(lon), lon_sat)
    valido = np.isfinite(lat) & np.isfinite(bcm) & (z < VZA_MAX)
    nube = valido & (bcm >= 1)
    cod = np.zeros(lat.shape, np.uint8)
    cod[valido & ~nube] = 1
    # "probablemente despejado" solo cuando la probabilidad de nube es dudosa (bruma, cirro tenue)
    # De noche NOAA etiqueta casi todo lo despejado como ACM=1, así que no se usa esa categoría.
    cod[valido & ~nube & (prob >= PROB_DUDOSA)] = 2
    h = np.where(np.isfinite(altura) & (altura > 0), altura, 4000.0)
    cod[nube & (h < H_BAJA)] = 3
    cod[nube & (h >= H_BAJA) & (h <= H_ALTA)] = 4
    cod[nube & (h > H_ALTA)] = 5

    n_niebla = 0
    if bt11 is not None and bt39 is not None:
        elev = elevacion_solar(np.nan_to_num(lat), np.nan_to_num(lon), t)
        btd = bt11 - bt39
        niebla = valido & ~nube & (elev < -5) & (btd > UMBRAL_NIEBLA_BTD) & (bt11 > BT11_MIN_NIEBLA)
        cod[niebla] = 6
        n_niebla = int(niebla.sum())

    # ── paralaje: desplazamiento de los píxeles nubosos ──
    dlat = np.zeros(lat.shape, np.float64)
    dlon = np.zeros(lat.shape, np.float64)
    if nube.any():
        a, b = desplazamiento_paralaje(lat[nube], lon[nube], h[nube], lon_sat)
        dlat[nube] = a
        dlon[nube] = b

    # ── rasterizado con submuestreo 3×3 por píxel ──
    g_r_lat = np.gradient(lat, axis=0)
    g_c_lat = np.gradient(lat, axis=1)
    g_r_lon = np.gradient(lon, axis=0)
    g_c_lon = np.gradient(lon, axis=1)
    rejilla = np.zeros((filas, cols), np.uint8)
    m = cod > 0
    base_lat, base_lon, c_m = lat[m] + dlat[m], lon[m] + dlon[m], cod[m]
    grl, gcl, grn, gcn = g_r_lat[m], g_c_lat[m], g_r_lon[m], g_c_lon[m]
    for fa, fb in itertools.product((-1 / 3, 0, 1 / 3), repeat=2):
        la = base_lat + fa * grl + fb * gcl
        lo = base_lon + fa * grn + fb * gcn
        r = np.floor((region["lat_max"] - la) / region["paso"]).astype(np.int64)
        c = np.floor((lo - region["lon_min"]) / region["paso"]).astype(np.int64)
        ok = (r >= 0) & (r < filas) & (c >= 0) & (c < cols) & np.isfinite(la) & np.isfinite(lo)
        np.maximum.at(rejilla, (r[ok], c[ok]), c_m[ok])

    # rellena huecos aislados con el vecino
    hueco = rejilla == 0
    if hueco.any():
        vec = np.zeros_like(rejilla)
        for dr, dc in ((1, 0), (-1, 0), (0, 1), (0, -1)):
            vec = np.maximum(vec, np.roll(np.roll(rejilla, dr, 0), dc, 1))
        rejilla[hueco] = vec[hueco]

    est = {k: int((rejilla == k).sum()) for k in CODIGOS}
    log(f"  celdas: " + ", ".join(f"{CODIGOS[k]}={v}" for k, v in est.items()) + f" · niebla px={n_niebla}")
    return rejilla, lon_sat


# ─────────────────────────── histórico y salida ───────────────────────────

def bajar_historico(sitio, dir_f):
    if not sitio:
        return []
    import requests
    base = sitio.rstrip("/") + "/data/"
    try:
        r = requests.get(base + "latest.json", timeout=20)
        r.raise_for_status()
        meta = r.json()
    except Exception as e:
        log(f"  (sin histórico previo: {e})")
        return []
    ok = []
    for c in meta.get("cuadros", []):
        destino = os.path.join(dir_f, os.path.basename(c["f"]))
        try:
            rr = requests.get(base + c["f"], timeout=30)
            rr.raise_for_status()
            with open(destino, "wb") as fh:
                fh.write(rr.content)
            ok.append(c)
        except Exception as e:
            log(f"  (no se pudo recuperar {c['f']}: {e})")
    log(f"  histórico recuperado: {len(ok)} cuadros")
    return ok


def guardar(rejilla, ruta):
    with open(ruta, "wb") as fh:
        fh.write(gzip.compress(rejilla.tobytes(), compresslevel=9, mtime=0))


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--salida", default="_site")
    ap.add_argument("--sitio", default=os.environ.get("SITIO_URL", ""),
                    help="URL pública del sitio ya desplegado (para recuperar el histórico)")
    ap.add_argument("--raiz", default="noaa-goes19")
    ap.add_argument("--local", action="store_true", help="--raiz es un directorio local (pruebas)")
    ap.add_argument("--cuadros", type=int, default=1, help="cuántas imágenes nuevas procesar (relleno inicial)")
    ap.add_argument("--paso-min", type=int, default=15)
    args = ap.parse_args(argv)

    dir_data = os.path.join(args.salida, "data")
    dir_f = os.path.join(dir_data, "f")
    os.makedirs(dir_f, exist_ok=True)

    cuadros = {c["t"]: c for c in bajar_historico(args.sitio, dir_f)}
    goes = GOES(args.raiz, local=args.local)
    ahora = datetime.now(timezone.utc)

    lista = goes.listar("ABI-L2-ACMF", ahora - timedelta(hours=HISTORICO_HORAS), ahora)
    if not lista:
        log("No hay máscaras ACMF recientes.")
    # si no hay histórico, rellena hacia atrás para que la tendencia tenga sentido
    n_nuevos = args.cuadros if cuadros else max(args.cuadros, 6)
    elegidos, ultimo = [], None
    for ts, f in reversed(lista):
        if ultimo is None or (ultimo - ts) >= timedelta(minutes=args.paso_min - 1):
            elegidos.append((ts, f))
            ultimo = ts
        if len(elegidos) >= n_nuevos:
            break

    lon_sat = -75.2
    for ts, f in reversed(elegidos):
        clave = ts.strftime("%Y-%m-%dT%H:%MZ")
        if clave in cuadros:
            log(f"{clave} ya procesado")
            continue
        log(f"Procesando {clave}")
        try:
            rejilla, lon_sat = procesar_cuadro(goes, ts, f)
        except Exception as e:
            log(f"  ERROR: {e}")
            continue
        nombre = f"f/{ts:%Y%m%d%H%M}.bin.gz"
        guardar(rejilla, os.path.join(dir_data, nombre))
        cuadros[clave] = {"t": clave, "f": nombre}

    # poda del histórico
    limite = (ahora - timedelta(hours=HISTORICO_HORAS)).strftime("%Y-%m-%dT%H:%MZ")
    orden = sorted((c for c in cuadros.values() if c["t"] >= limite), key=lambda c: c["t"])[-MAX_CUADROS:]
    vivos = {os.path.basename(c["f"]) for c in orden}
    for n in os.listdir(dir_f):
        if n not in vivos:
            os.remove(os.path.join(dir_f, n))

    filas = int(round((REGION["lat_max"] - REGION["lat_min"]) / REGION["paso"]))
    cols = int(round((REGION["lon_max"] - REGION["lon_min"]) / REGION["paso"]))
    try:
        from movimiento import campo_movimiento
        mov = campo_movimiento(orden, dir_data, filas, cols, REGION, log=log)
    except Exception as e:
        log(f"  movimiento: ERROR {e}")
        mov = None
    meta = {
        "generado": ahora.strftime("%Y-%m-%dT%H:%M:%SZ"),
        "satelite": "GOES-19",
        "lon_sat": lon_sat,
        "rejilla": {**REGION, "filas": filas, "cols": cols},
        "codigos": CODIGOS,
        "cuadros": orden,
        "movimiento": mov,
    }
    with open(os.path.join(dir_data, "latest.json"), "w", encoding="utf-8") as fh:
        json.dump(meta, fh, ensure_ascii=False, indent=1)
    log(f"Listo: {len(orden)} cuadros en el histórico")
    return 0 if orden else 1


if __name__ == "__main__":
    sys.exit(main())
