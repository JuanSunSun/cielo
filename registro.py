"""Registro de pronósticos y resultados para entrenar la corrección estadística.

Cada hora (como mucho) guarda, para ~500 puntos fijos, lo que el pronóstico por advección sabe en ese
momento (nubosidad actual, nubosidad extrapolada a +1/+2/+3 h, calidad del movimiento) y, cuando pasan
1, 2 y 3 horas, rellena lo que de verdad ocurrió en ese punto. Con eso se entrena `entrenar.py`.

Formato: un archivo por día, `registro/AAAA-MM-DD.jsonl.gz`; cada línea es una pasada con columnas
paralelas (un valor por punto, null si no hay dato). Tope de tamaño total con rotación de lo más antiguo.

Uso:
    python registro.py --sitio https://juansunsun.github.io/cielo --datos datos
    python registro.py --local DIR_CON_latest.json_Y_f --datos datos      (pruebas)
"""
import argparse
import gzip
import json
import math
import os
import sys
import tempfile
from datetime import datetime, timedelta, timezone

import numpy as np

PASO_GRADOS = 1.0            # separación de la malla de puntos
RADIO_KM = 10.0              # radio con el que se mide la nubosidad
HORIZONTES_H = (1, 2, 3)
INTERVALO_MIN = 50           # no se registra más de una vez cada ~hora
TOL_ETIQUETA_MIN = 14        # un cuadro vale como «observado» si está a menos de esto del instante objetivo
MAX_BYTES = 40 * 1024 * 1024 # tope del directorio registro/ (se borra lo más antiguo)
NUBE = (3, 4, 5, 6)


def log(*a):
    print(*a, file=sys.stderr, flush=True)


# ───────────────────────────── utilidades ─────────────────────────────

def t_de(s):
    return datetime.strptime(s, "%Y-%m-%dT%H:%MZ").replace(tzinfo=timezone.utc)


def s_de(t):
    return t.strftime("%Y-%m-%dT%H:%MZ")


def elev_solar(lat, lon, t):
    """Elevación solar en grados (misma fórmula que la app)."""
    doy = t.timetuple().tm_yday
    h = t.hour + t.minute / 60
    g = 2 * math.pi / 365 * (doy - 1 + (h - 12) / 24)
    dec = (.006918 - .399912 * math.cos(g) + .070257 * math.sin(g) - .006758 * math.cos(2 * g)
           + .000907 * math.sin(2 * g) - .002697 * math.cos(3 * g) + .00148 * math.sin(3 * g))
    eq = 229.18 * (.000075 + .001868 * math.cos(g) - .032077 * math.sin(g) - .014615 * math.cos(2 * g)
                   - .040849 * math.sin(2 * g))
    ha = math.radians(((h * 60 + eq + 4 * lon) / 4 - 180))
    la = math.radians(lat)
    c = math.sin(la) * math.sin(dec) + math.cos(la) * math.cos(dec) * math.cos(ha)
    return 90 - math.degrees(math.acos(max(-1, min(1, c))))


def puntos(G):
    """Malla fija de puntos dentro de la rejilla (lat, lon)."""
    out = []
    lat = G["lat_max"] - PASO_GRADOS / 2
    while lat > G["lat_min"] + PASO_GRADOS / 2 - 1e-9:
        lon = G["lon_min"] + PASO_GRADOS / 2
        while lon < G["lon_max"] - PASO_GRADOS / 2 + 1e-9:
            out.append((round(lat, 3), round(lon, 3)))
            lon += PASO_GRADOS
        lat -= PASO_GRADOS
    return out


def fraccion_disco(g, G, lat, lon, radio, busq=0):
    """(fracción de nube, píxeles válidos, distancia a la nube más cercana en `busq` km o None).
    Misma geometría que `analizar()` de la app; g es uint8 (filas, cols)."""
    paso = G["paso"]
    kx = 111.32 * math.cos(math.radians(lat)) * paso
    ky = 110.57 * paso
    r0 = (G["lat_max"] - lat) / paso - .5
    c0 = (lon - G["lon_min"]) / paso - .5
    if r0 < 0 or r0 >= G["filas"] or c0 < 0 or c0 >= G["cols"]:
        return None, 0, None
    R = max(radio, busq)
    dr, dc = int(math.ceil(R / ky)), int(math.ceil(R / kx))
    ra, rb = max(0, int(math.floor(r0 - dr))), min(G["filas"] - 1, int(math.ceil(r0 + dr)))
    ca, cb = max(0, int(math.floor(c0 - dc))), min(G["cols"] - 1, int(math.ceil(c0 + dc)))
    sub = g[ra:rb + 1, ca:cb + 1]
    rr, cc = np.mgrid[ra:rb + 1, ca:cb + 1]
    d = np.hypot((cc - c0) * kx, (r0 - rr) * ky)
    dentro = d <= radio
    valido = dentro & (sub > 0)
    n = int(valido.sum())
    if n == 0:
        v = g[int(round(r0)), int(round(c0))]
        frac, n = (None if v == 0 else float(v >= 3)), int(v > 0)
    else:
        frac = float((valido & (sub >= 3)).sum() / n)
    cerca = None
    if busq:
        nube = (sub >= 3) & (d <= busq)
        if nube.any():
            cerca = float(d[nube].min())
    return frac, n, cerca


def vector_local(M, lat, lon):
    """Mediana ponderada de los vectores a < 450 km (misma lógica que la app).
    Devuelve (u, v, dispersión, coherencia) o None."""
    if not M or not M.get("vectores"):
        return None
    kx = 111.32 * math.cos(math.radians(lat))
    cand = []
    for x in M["vectores"]:
        d = math.hypot((x["lat"] - lat) * 110.57, (x["lon"] - lon) * kx)
        if d <= 450:
            cand.append((x["u"], x["v"], 1 if x.get("k") else 0, x["c"] / (1 + (d / 200) ** 2)))
    if len(cand) < 2:
        g = M.get("global")
        return (g["u"], g["v"], None, None) if g else None

    def mediana(vals, w):
        o = np.argsort(vals)
        cw = np.cumsum(np.array(w)[o])
        return float(np.array(vals)[o][np.searchsorted(cw, cw[-1] / 2)])

    w = [c[3] for c in cand]
    u, v = mediana([c[0] for c in cand], w), mediana([c[1] for c in cand], w)
    disp = mediana([math.hypot(c[0] - u, c[1] - v) for c in cand], w)
    coh = sum(c[2] * c[3] for c in cand) / (sum(w) or 1)
    return u, v, disp, coh


def redondea(x, n=3):
    return None if x is None or (isinstance(x, float) and math.isnan(x)) else round(float(x), n)


# ─────────────────────────── características ───────────────────────────

def caracteristicas(meta, cuadros, pts):
    """Columnas por punto para el cuadro más reciente. cuadros: {clave_t: array}."""
    G = meta["rejilla"]
    ultimo = meta["cuadros"][-1]
    g = cuadros[ultimo["t"]]
    M = meta.get("movimiento")
    cols = {k: [] for k in ("f0", "a1", "a2", "a3", "u", "v", "d", "k", "c")}
    for lat, lon in pts:
        f0, n, cerca = fraccion_disco(g, G, lat, lon, RADIO_KM, busq=150)
        cols["f0"].append(redondea(f0, 3)); cols["c"].append(redondea(cerca, 1))
        vec = vector_local(M, lat, lon) if f0 is not None else None
        if vec is None:
            for k in ("a1", "a2", "a3", "u", "v", "d", "k"): cols[k].append(None)
            continue
        u, v, disp, coh = vec
        cols["u"].append(redondea(u, 1)); cols["v"].append(redondea(v, 1))
        cols["d"].append(redondea(disp, 1)); cols["k"].append(redondea(coh, 2))
        kx = 111.32 * math.cos(math.radians(lat))
        vel = math.hypot(u, v)
        for h in HORIZONTES_H:
            la, lo = lat - v * h / 110.57, lon - u * h / kx
            r_ef = math.hypot(RADIO_KM, 0.2 * vel * h)
            fa, _, _ = fraccion_disco(g, G, la, lo, r_ef)
            cols[f"a{h}"].append(redondea(fa, 3))
    return cols


def observado(meta, cuadros, pts, t_objetivo):
    """Fracción observada en cada punto en el cuadro más cercano a t_objetivo (o None si no hay uno cerca)."""
    G = meta["rejilla"]
    mejor = None
    for c in meta["cuadros"]:
        dif = abs((t_de(c["t"]) - t_objetivo).total_seconds()) / 60
        if mejor is None or dif < mejor[0]:
            mejor = (dif, c["t"])
    if mejor is None or mejor[0] > TOL_ETIQUETA_MIN:
        return None
    g = cuadros[mejor[1]]
    return [redondea(fraccion_disco(g, G, lat, lon, RADIO_KM)[0], 3) for lat, lon in pts]


# ───────────────────────────── almacenamiento ─────────────────────────────

def ruta_dia(dir_reg, t):
    return os.path.join(dir_reg, f"{t:%Y-%m-%d}.jsonl.gz")


def leer_dia(ruta):
    if not os.path.exists(ruta):
        return []
    with gzip.open(ruta, "rt", encoding="utf-8") as fh:
        return [json.loads(l) for l in fh if l.strip()]


def escribir_dia(ruta, filas):
    tmp = ruta + ".tmp"
    with gzip.open(tmp, "wt", encoding="utf-8", compresslevel=9) as fh:
        for f in filas:
            fh.write(json.dumps(f, separators=(",", ":")) + "\n")
    os.replace(tmp, ruta)


def rotar(dir_reg, max_bytes=MAX_BYTES):
    """Borra los días más antiguos hasta que el directorio quepa en max_bytes."""
    arch = sorted(f for f in os.listdir(dir_reg) if f.endswith(".jsonl.gz"))
    tam = lambda: sum(os.path.getsize(os.path.join(dir_reg, f)) for f in arch)
    borrados = []
    while len(arch) > 1 and tam() > max_bytes:
        os.remove(os.path.join(dir_reg, arch[0])); borrados.append(arch.pop(0))
    return borrados, tam()


def registrar(meta, cuadros, dir_datos, ahora=None):
    """Añade una pasada y rellena las etiquetas pendientes. Devuelve un resumen (dict)."""
    dir_reg = os.path.join(dir_datos, "registro")
    os.makedirs(dir_reg, exist_ok=True)
    G = meta["rejilla"]
    pts = puntos(G)
    ruta_pts = os.path.join(dir_reg, "puntos.json")
    if os.path.exists(ruta_pts):
        guardados = [tuple(p) for p in json.load(open(ruta_pts))]
        if guardados != pts:
            raise SystemExit("La malla de puntos cambió: cambia PASO_GRADOS con cuidado o empieza un registro nuevo")
    else:
        json.dump(pts, open(ruta_pts, "w"))

    t0 = t_de(meta["cuadros"][-1]["t"])
    estado_ruta = os.path.join(dir_datos, "estado.json")
    estado = json.load(open(estado_ruta)) if os.path.exists(estado_ruta) else {}
    resumen = {"t0": s_de(t0), "etiquetas": 0, "nueva": False}

    # 1) etiquetas pendientes: pasadas de hace 1, 2 y 3 h (también del día anterior)
    dias = {}
    for dd in (t0, t0 - timedelta(days=1)):
        r = ruta_dia(dir_reg, dd)
        dias[r] = leer_dia(r)
    tocados = set()
    for r, filas in dias.items():
        for fila in filas:
            tf = t_de(fila["t"])
            for h in HORIZONTES_H:
                clave = f"o{h}"
                if fila.get(clave) is not None:
                    continue
                objetivo = tf + timedelta(hours=h)
                if objetivo > t0 + timedelta(minutes=TOL_ETIQUETA_MIN):
                    continue
                obs = observado(meta, cuadros, pts, objetivo)
                if obs is not None:
                    fila[clave] = obs; tocados.add(r); resumen["etiquetas"] += 1

    # 2) pasada nueva (como mucho una por hora)
    ultimo = estado.get("ultimo")
    if ultimo is None or (t0 - t_de(ultimo)).total_seconds() / 60 >= INTERVALO_MIN:
        cols = caracteristicas(meta, cuadros, pts)
        fila = {"t": s_de(t0), **cols, "o1": None, "o2": None, "o3": None}
        r = ruta_dia(dir_reg, t0)
        dias.setdefault(r, leer_dia(r)).append(fila)
        tocados.add(r)
        estado["ultimo"] = s_de(t0)
        resumen["nueva"] = True

    for r in tocados:
        escribir_dia(r, dias[r])
    borrados, tam = rotar(dir_reg)
    estado["bytes"] = tam
    json.dump(estado, open(estado_ruta, "w"))
    resumen.update(bytes=tam, borrados=borrados, puntos=len(pts))
    return resumen


# ──────────────────────────────── entrada ────────────────────────────────

def cargar_cuadros(meta, leer):
    out = {}
    for c in meta["cuadros"]:
        raw = gzip.decompress(leer(c["f"]))
        G = meta["rejilla"]
        out[c["t"]] = np.frombuffer(raw, np.uint8).reshape(G["filas"], G["cols"])
    return out


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--sitio", default="")
    ap.add_argument("--local", default="", help="directorio con latest.json y f/ (pruebas)")
    ap.add_argument("--datos", default="datos")
    a = ap.parse_args(argv)
    if a.local:
        meta = json.load(open(os.path.join(a.local, "latest.json")))
        cuadros = cargar_cuadros(meta, lambda f: open(os.path.join(a.local, f), "rb").read())
    else:
        import requests
        base = a.sitio.rstrip("/") + "/data/"
        meta = requests.get(base + "latest.json", timeout=30).json()
        cuadros = cargar_cuadros(meta, lambda f: requests.get(base + f, timeout=60).content)
    r = registrar(meta, cuadros, a.datos)
    log(json.dumps(r, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
