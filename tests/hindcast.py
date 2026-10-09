"""Verificación retrospectiva del movimiento: ¿predecir con los vectores supera a «no cambia nada»?

Uso: python tests/hindcast.py DIRECTORIO_CON_latest.json_Y_f/
Para cada cuadro T[i] de una cadena a 30 min, calcula los vectores con (T[i-1], T[i]) y
predice T[i+h] (h = 30, 60, 90 min) moviendo T[i]. Compara con el cuadro real y con la persistencia.
"""
import json, math, os, sys
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import movimiento as M

d = sys.argv[1]
meta = json.load(open(os.path.join(d, "latest.json")))
G = meta["rejilla"]; F, C, PASO_G = G["filas"], G["cols"], G["paso"]
cuadros = {c["t"]: M._leer(os.path.join(d, c["f"]), F, C) for c in meta["cuadros"]}
ts = sorted(cuadros)
from datetime import datetime, timedelta
T = lambda s: datetime.strptime(s, "%Y-%m-%dT%H:%MZ")
S = lambda t: t.strftime("%Y-%m-%dT%H:%MZ")

def campo_denso(vs, a_shape):
    nr = (F - M.TESELA) // M.PASO + 1; nc = (C - M.TESELA) // M.PASO + 1
    U = np.zeros((nr, nc), np.float32); V = np.zeros_like(U); W = np.zeros_like(U)
    for x in vs:
        r = int(round(((G["lat_max"] - x["lat"]) / PASO_G - M.TESELA / 2) / M.PASO))
        c = int(round(((x["lon"] - G["lon_min"]) / PASO_G - M.TESELA / 2) / M.PASO))
        if 0 <= r < nr and 0 <= c < nc: U[r, c], V[r, c], W[r, c] = x["u"], x["v"], 1
    # bilinear sobre centros de tesela; las teselas sin vector pesan 0 y se rellenan con vecinos (suavizado normalizado)
    def interp(A):
        rr = (np.arange(F) - M.TESELA / 2) / M.PASO; cc = (np.arange(C) - M.TESELA / 2) / M.PASO
        r0 = np.clip(np.floor(rr).astype(int), 0, nr - 2); c0 = np.clip(np.floor(cc).astype(int), 0, nc - 2)
        fr = np.clip(rr - r0, 0, 1)[:, None]; fc = np.clip(cc - c0, 0, 1)[None, :]
        g = lambda M_: (M_[r0][:, c0] * (1 - fr) * (1 - fc) + M_[r0 + 1][:, c0] * fr * (1 - fc)
                        + M_[r0][:, c0 + 1] * (1 - fr) * fc + M_[r0 + 1][:, c0 + 1] * fr * fc)
        return g(A * W), g(W)
    nu, w = interp(U); nv, _ = interp(V)
    w = np.maximum(w, 1e-6)
    return nu / w, nv / w, interp(W)[0]    # u, v (km/h) y peso

lat = G["lat_max"] - (np.arange(F) + .5) * PASO_G
kx_px = (111.32 * np.cos(np.radians(lat)) * PASO_G)[:, None]; ky_px = 110.57 * PASO_G
rows, cols = np.mgrid[0:F, 0:C]

def avanzar(frame, u, v, w, horas):
    drow = -v * horas / ky_px; dcol = u * horas / kx_px
    sr = np.clip(np.rint(rows - drow).astype(int), 0, F - 1); sc = np.clip(np.rint(cols - dcol).astype(int), 0, C - 1)
    out = frame[sr, sc]
    # donde no hay información de movimiento (peso ~ 0) se deja la persistencia
    return np.where(w > 0.15, out, frame)

def caja(m, k=9):
    """Media móvil k×k (fracción de nube en ~20 km con celdas de 2 km)."""
    m = m.astype(np.float32); pad = k // 2
    c = np.cumsum(np.cumsum(np.pad(m, pad, mode="edge"), 0), 1)
    c = np.pad(c, ((1, 0), (1, 0)))
    return (c[k:, k:] - c[:-k, k:] - c[k:, :-k] + c[:-k, :-k]) / (k * k)

clase = lambda f: np.digitize(f, [0.10, 0.50])      # 0 despejado, 1 parcial, 2 nublado
cn = lambda x: x >= 3
SUB = (slice(8, F, 14), slice(8, C, 14))
res = {}
for t0 in ts:
    tprev = S(T(t0) - timedelta(minutes=30))
    if tprev not in cuadros: continue
    a, b = cuadros[tprev], cuadros[t0]
    vs = M.campo_desde_cuadros(a, b, 30, {"lat_max": G["lat_max"], "lon_min": G["lon_min"], "paso": PASO_G})
    u, v, w = campo_denso(vs, a.shape)
    for h in (0.5, 1.0, 1.5):
        th = S(T(t0) + timedelta(minutes=int(h * 60)))
        if th not in cuadros: continue
        base, real = cuadros[t0], cuadros[th]
        ok = ((real > 0) & (base > 0))[SUB]
        f_b = caja(cn(base))[SUB]; f_r = caja(cn(real))[SUB]
        f_p = caja(cn(avanzar(base, u, v, w, h)))[SUB]
        cb, cr, cp = clase(f_b)[ok], clase(f_r)[ok], clase(f_p)[ok]
        cambia = cb != cr
        r = dict(mae_pers=np.abs(f_b - f_r)[ok].mean(), mae_adv=np.abs(f_p - f_r)[ok].mean(),
                 acc_pers=(cb == cr).mean(), acc_adv=(cp == cr).mean(),
                 n_cambia=cambia.mean(),
                 hit=((cp == cr) & cambia).sum() / max(cambia.sum(), 1),            # cambios reales bien predichos
                 falsa=((cp != cb) & ~cambia).sum() / max((~cambia).sum(), 1))      # cambios predichos que no ocurrieron
        res.setdefault(h, []).append(r)
print("Fracción de nubes en ~20 km (la magnitud que usa la app), 3 clases: despejado ≤10 %, parcial, nublado >50 %\n")
for h, L in sorted(res.items()):
    m = {k: np.mean([x[k] for x in L]) for k in L[0]}
    print(f"+{int(h*60):3d} min (n={len(L)}):  error medio de la fracción  persistencia {m['mae_pers']*100:5.2f} → advección {m['mae_adv']*100:5.2f}  "
          f"({(1-m['mae_adv']/m['mae_pers'])*100:3.0f} % mejor)\n"
          f"            clase correcta  persistencia {m['acc_pers']*100:5.1f} % → advección {m['acc_adv']*100:5.1f} %\n"
          f"            hay cambio de clase en el {m['n_cambia']*100:4.1f} % de los puntos; la advección acierta {m['hit']*100:4.1f} % de ellos "
          f"y da falsa alarma de cambio en {m['falsa']*100:4.1f} %")
