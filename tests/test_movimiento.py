"""Verifica movimiento.py con un campo de nubes sintético desplazado una cantidad conocida."""
import math, sys, os
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import movimiento as M

REG = dict(lat_max=-17.5, lat_min=-56.0, lon_min=-76.0, lon_max=-64.0, paso=0.02)
F, C = 1925, 600

def campo(seed):
    rng = np.random.default_rng(seed)
    x = rng.standard_normal((F + 200, C + 200)).astype(np.float32)
    fy = np.fft.fftfreq(x.shape[0])[:, None]; fx = np.fft.fftfreq(x.shape[1])[None, :]
    k = np.exp(-(fx ** 2 + fy ** 2) / (2 * 0.012 ** 2))      # estructuras de ~15-25 px
    s = np.fft.ifft2(np.fft.fft2(x) * k).real
    return s

def cuadro(s, dy, dx):
    """Recorta s desplazado (dy filas hacia abajo, dx columnas a la derecha)."""
    ys, xs = 100 - int(round(dy)), 100 - int(round(dx))
    sub = s[ys:ys + F, xs:xs + C]
    u = (sub > np.percentile(sub, 55)).astype(np.uint8)
    cod = np.where(u == 1, 4, 1).astype(np.uint8)
    return cod

def prueba(nombre, dx_px, dy_px, dt_min=30):
    s = campo(1)
    a = cuadro(s, 0, 0); b = cuadro(s, dy_px, dx_px)
    vs = M.campo_desde_cuadros(a, b, dt_min, REG)
    ky = 110.57 * REG["paso"]
    esp_v = -dy_px * ky / (dt_min / 60)
    errs = []
    for x in vs:
        kx = 111.32 * math.cos(math.radians(x["lat"])) * REG["paso"]
        esp_u = dx_px * kx / (dt_min / 60)
        errs.append(math.hypot(x["u"] - esp_u, x["v"] - esp_v))
    g = M.resumen_global(vs)
    print(f"{nombre}: {len(vs)} vectores; error mediano {np.median(errs):.1f} km/h, p90 {np.percentile(errs, 90):.1f}; "
          f"esperado v={esp_v:.1f}; global={g}")
    return len(vs), float(np.median(errs))

ok = True
for nombre, dx, dy in [("hacia el este 9 px", 9, 0), ("hacia el SE (NO→SE)", 6, 5), ("hacia el norte", 0, -8), ("lento 2 px", 2, 1)]:
    n, e = prueba(nombre, dx, dy)
    if n < 100 or e > 3.0:
        ok = False
print("OK" if ok else "FALLO")
sys.exit(0 if ok else 1)
