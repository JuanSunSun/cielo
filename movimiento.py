"""Campo de movimiento de las nubes a partir de dos cuadros de la máscara.

1. Correlación de fase por teselas (ventana de Hann, filtro paso bajo, pico subpíxel).
2. Problema de la apertura: un borde recto de nubes solo fija la componente perpendicular a
   él. Esas teselas aportan una restricción (la componente normal) en lugar de un vector.
3. Se resuelve un vector por tesela con mínimos cuadrados ponderados y reponderación robusta
   usando las restricciones de las teselas vecinas (< ~300 km).  Si las orientaciones de los
   bordes no bastan para fijar el vector (todos paralelos) esa zona se deja sin dato.

Solo usa numpy.

Convención: un cuadro es uint8 (filas N→S, columnas O→E) con códigos de cielo; nube = código >= 3.
u = este, v = norte, en km/h.  Una nube que estaba en x en el cuadro antiguo está en x + (u, v)·Δt
en el reciente.
"""
import gzip
import math
import os
from datetime import datetime

import numpy as np

TESELA = 96          # px (~210 km con celdas de 0,02°)
PASO = 48            # separación entre teselas
SIGMA_F = 0.07       # ciclos/px del filtro gaussiano paso bajo (~2,3 px de suavizado)
MAX_DESPL = 32       # px máximos buscados (≈ 140 km/h a 30 min)
FRAC_MIN, FRAC_MAX = 0.04, 0.96   # fracción de nube que da textura a la tesela
VALIDO_MIN = 0.85    # fracción mínima de píxeles con dato en ambos cuadros
CONF_MIN = 0.12      # altura mínima del pico de correlación
AISLAMIENTO_PX = 4   # distancia mínima al segundo máximo
RATIO_CRESTA = 0.80  # segundo/primer máximo a partir del cual el pico es una cresta (apertura)
VEL_MAX = 140.0      # km/h
RADIO_VECINOS_KM = 320.0
ESCALA_KM = 200.0    # decaimiento gaussiano del peso con la distancia
DT_IDEAL_MIN = 30
DT_MIN, DT_MAX = 14, 65


def _t(c):
    return datetime.strptime(c["t"], "%Y-%m-%dT%H:%MZ")


def _leer(ruta, filas, cols):
    with open(ruta, "rb") as fh:
        b = gzip.decompress(fh.read())
    return np.frombuffer(b, np.uint8).reshape(filas, cols)


def _ventana(n):
    h = np.hanning(n).astype(np.float32)
    return np.outer(h, h)


def _filtro(n):
    f = np.fft.fftfreq(n).astype(np.float32)
    fr = np.fft.rfftfreq(n).astype(np.float32)
    return np.exp(-(f[:, None] ** 2 + fr[None, :] ** 2) / (2 * SIGMA_F ** 2))


def _desplazamiento(a, b, win, filt):
    """(dx, dy, pico, ambiguo): desplazamiento de b respecto a a, en píxeles (dy hacia abajo)."""
    A = np.fft.rfft2((a - a.mean()) * win)
    B = np.fft.rfft2((b - b.mean()) * win)
    R = B * np.conj(A)
    R /= np.abs(R) + 1e-6
    R *= filt
    c = np.fft.irfft2(R, s=a.shape)
    n = a.shape[0]
    c = c / max(float(2 * filt.sum()) / (n * n), 1e-9)     # un desplazamiento perfecto da ≈ 1
    vi = np.r_[0:MAX_DESPL + 1, n - MAX_DESPL:n]
    sg = np.where(vi <= n // 2, vi, vi - n).astype(np.float32)
    sub = c[np.ix_(vi, vi)]
    jy, jx = np.unravel_index(int(np.argmax(sub)), sub.shape)
    iy, ix = int(vi[jy]), int(vi[jx])
    pico = float(c[iy, ix])
    d2 = (sg[:, None] - sg[jy]) ** 2 + (sg[None, :] - sg[jx]) ** 2
    lejos = sub[d2 > AISLAMIENTO_PX ** 2]
    segundo = float(lejos.max()) if lejos.size else 0.0
    ambiguo = bool(pico > 0 and segundo / pico > RATIO_CRESTA)

    def parab(m, o, p):
        den = m - 2 * o + p
        return 0.0 if abs(den) < 1e-12 else 0.5 * (m - p) / den

    sy = parab(c[(iy - 1) % n, ix], c[iy, ix], c[(iy + 1) % n, ix])
    sx = parab(c[iy, (ix - 1) % n], c[iy, ix], c[iy, (ix + 1) % n])
    dy = (iy + sy) if iy <= n // 2 else (iy + sy - n)
    dx = (ix + sx) if ix <= n // 2 else (ix + sx - n)
    return float(dx), float(dy), pico, ambiguo


def _normal(t, kx, ky):
    """Dirección (este, norte) del gradiente dominante de la tesela, unitaria."""
    gy, gx = np.gradient(t)
    ge, gn = gx / kx, -gy / ky                 # gradiente por km (este, norte)
    sxx, snn, sxn = float((ge * ge).sum()), float((gn * gn).sum()), float((ge * gn).sum())
    ang = 0.5 * math.atan2(2 * sxn, sxx - snn)
    return math.cos(ang), math.sin(ang)


def elegir_pareja(cuadros):
    """El cuadro más reciente y el más cercano a 30 min antes (entre 14 y 65 min)."""
    if len(cuadros) < 2:
        return None
    nuevo = cuadros[-1]
    mejor = None
    for c in cuadros[:-1]:
        dt = (_t(nuevo) - _t(c)).total_seconds() / 60
        if DT_MIN <= dt <= DT_MAX:
            d = abs(dt - DT_IDEAL_MIN)
            if mejor is None or d < mejor[0]:
                mejor = (d, c, dt)
    return None if mejor is None else (mejor[1], nuevo, mejor[2])


def _restricciones(a, b, dt_min, region):
    """Una entrada por tesela: posición y restricciones lineales n·(u,v) = s con su peso."""
    filas, cols = a.shape
    paso = region["paso"]
    ma, mb = (a >= 3).astype(np.float32), (b >= 3).astype(np.float32)
    va, vb = (a > 0), (b > 0)
    win, filt = _ventana(TESELA), _filtro(TESELA)
    ky = 110.57 * paso
    horas = dt_min / 60.0
    out = []
    for r in range(0, filas - TESELA + 1, PASO):
        for c in range(0, cols - TESELA + 1, PASO):
            sl = (slice(r, r + TESELA), slice(c, c + TESELA))
            if (va[sl] & vb[sl]).mean() < VALIDO_MIN:
                continue
            ta, tb = ma[sl], mb[sl]
            if not (FRAC_MIN < 0.5 * (ta.mean() + tb.mean()) < FRAC_MAX):
                continue
            dx, dy, pico, amb = _desplazamiento(ta, tb, win, filt)
            if pico < CONF_MIN:
                continue
            lat = region["lat_max"] - (r + TESELA / 2) * paso
            lon = region["lon_min"] + (c + TESELA / 2) * paso
            kx = 111.32 * math.cos(math.radians(lat)) * paso
            u, v = dx * kx / horas, -dy * ky / horas
            if math.hypot(u, v) > VEL_MAX * 1.5:
                continue
            w = min(1.0, pico)
            if amb:
                ne, nn = _normal(0.5 * (ta + tb), kx, ky)
                cons = [(ne, nn, u * ne + v * nn, 0.6 * w)]
            else:
                cons = [(1.0, 0.0, u, w), (0.0, 1.0, v, w)]
            out.append({"lat": lat, "lon": lon, "cons": cons})
    return out


def _resolver(tesela, todas):
    """(u, v, calidad) en una tesela con las restricciones vecinas, o None si no está determinado."""
    kx = 111.32 * math.cos(math.radians(tesela["lat"]))
    N, S, W = [], [], []
    for o in todas:
        d = math.hypot((o["lat"] - tesela["lat"]) * 110.57, (o["lon"] - tesela["lon"]) * kx)
        if d > RADIO_VECINOS_KM:
            continue
        pd = math.exp(-(d / ESCALA_KM) ** 2)
        for ne, nn, s, w in o["cons"]:
            N.append((ne, nn)); S.append(s); W.append(w * pd)
    if len(N) < 3:
        return None
    N, S, W = np.array(N), np.array(S), np.array(W)
    w = W.copy()
    uv, ev = None, None
    for _ in range(4):      # reponderación robusta (Cauchy, 8 km/h)
        A = (N * w[:, None]).T @ N
        ev = np.linalg.eigvalsh(A)
        if ev[0] < 0.12 * max(w.sum(), 1e-9) or w.sum() < 0.8:
            return None
        uv = np.linalg.solve(A, (N * w[:, None]).T @ S)
        w = W / (1.0 + ((N @ uv - S) / 8.0) ** 2)
    calidad = float(min(1.0, w.sum() / 4.0) * min(1.0, ev[0] / (0.35 * w.sum())))
    return float(uv[0]), float(uv[1]), calidad


def campo_desde_cuadros(a, b, dt_min, region):
    """a: cuadro antiguo, b: reciente (uint8). Devuelve la lista de vectores."""
    todas = _restricciones(a, b, dt_min, region)
    out = []
    for t in todas:
        r = _resolver(t, todas)
        if r is None or math.hypot(r[0], r[1]) > VEL_MAX or r[2] < 0.15:
            continue
        out.append({"lat": round(t["lat"], 3), "lon": round(t["lon"], 3),
                    "u": round(r[0], 1), "v": round(r[1], 1), "c": round(r[2], 2)})
    return out


def resumen_global(vs):
    if not vs:
        return None
    w = np.array([x["c"] for x in vs])
    out = {}
    for k in ("u", "v"):
        a = np.array([x[k] for x in vs])
        o = np.argsort(a)
        cw = np.cumsum(w[o])
        out[k] = round(float(a[o][np.searchsorted(cw, cw[-1] / 2)]), 1)
    out["n"] = len(vs)
    return out


def _marcar_coherencia(vs, vs_previo, tol=8.0):
    """Añade k=1 a los vectores que coinciden (< tol km/h) con el medido en la pareja anterior."""
    previo = {(x["lat"], x["lon"]): x for x in vs_previo}
    for x in vs:
        o = previo.get((x["lat"], x["lon"]))
        x["k"] = 1 if o is not None and math.hypot(x["u"] - o["u"], x["v"] - o["v"]) < tol else 0
    return vs


def campo_movimiento(orden, dir_data, filas, cols, region, log=print):
    """orden: lista de cuadros {t,f} ordenada; devuelve el dict para latest.json o None."""
    pareja = elegir_pareja(orden)
    if not pareja:
        log("  movimiento: no hay pareja de cuadros adecuada")
        return None
    ca, cb, dt = pareja
    a = _leer(os.path.join(dir_data, ca["f"]), filas, cols)
    b = _leer(os.path.join(dir_data, cb["f"]), filas, cols)
    vs = campo_desde_cuadros(a, b, dt, region)
    # coherencia temporal: ¿el movimiento de la pareja anterior (ca respecto a ~30 min antes) coincide?
    previa = elegir_pareja(orden[:orden.index(ca) + 1]) if ca in orden else None
    if previa:
        cc, _, dt2 = previa
        c = _leer(os.path.join(dir_data, cc["f"]), filas, cols)
        _marcar_coherencia(vs, campo_desde_cuadros(c, a, dt2, region))
    else:
        for x in vs:
            x["k"] = 0
    glob = resumen_global(vs)
    log(f"  movimiento: {len(vs)} vectores entre {ca['t']} y {cb['t']} (Δt {dt:.0f} min), "
        f"{sum(x['k'] for x in vs)} coherentes con la pareja anterior"
        + (f", mediana u={glob['u']} v={glob['v']} km/h" if glob else ""))
    return {"t0": ca["t"], "t1": cb["t"], "dt_min": round(dt, 1), "vectores": vs, "global": glob}
