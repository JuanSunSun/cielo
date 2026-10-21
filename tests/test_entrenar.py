"""Prueba del entrenador con un registro sintético donde se sabe que corregir ayuda
(la nube real se parece a 0,6·a_h + 0,3·f0 + ruido, y a_h sobrestima cuando la velocidad es baja)."""
import json, os, sys, tempfile, gzip
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import registro as R, entrenar as E

rng = np.random.default_rng(3)
G = dict(lat_max=-17.5, lat_min=-56.0, lon_min=-76.0, lon_max=-64.0)
pts = R.puntos(G); n = len(pts)
tmp = tempfile.mkdtemp(); reg = os.path.join(tmp, "registro"); os.makedirs(reg)
json.dump(pts, open(os.path.join(reg, "puntos.json"), "w"))
from datetime import datetime, timedelta, timezone
t = datetime(2026, 10, 1, 0, 20, tzinfo=timezone.utc)
dias = {}
for i in range(24 * 7):
    f0 = np.clip(rng.beta(.6, .6, n), 0, 1); vel = rng.gamma(2, 8, n)
    fila = {"t": R.s_de(t), "f0": f0.round(3).tolist(), "u": (vel * .8).round(1).tolist(), "v": (vel * .2).round(1).tolist(),
            "d": rng.uniform(0, 20, n).round(1).tolist(), "k": rng.uniform(0, 1, n).round(2).tolist(),
            "c": rng.uniform(0, 150, n).round(1).tolist()}
    for h in (1, 2, 3):
        a = np.clip(f0 + rng.normal(0, .25, n) * (h / 2), 0, 1)
        o = np.clip(.6 * a + .3 * f0 + rng.normal(0, .12, n) + .1 * (vel < 8), 0, 1)
        fila[f"a{h}"] = a.round(3).tolist(); fila[f"o{h}"] = o.round(3).tolist()
    dias.setdefault(f"{t:%Y-%m-%d}", []).append(fila)
    t += timedelta(hours=1)
for d, filas in dias.items():
    R.escribir_dia(os.path.join(reg, d + ".jsonl.gz"), filas)
sys.argv = ["entrenar.py", "--registro", reg, "--salida", os.path.join(tmp, "out")]
E.main()
m = json.load(open(os.path.join(tmp, "out", "modelo.json")))
ok = "1" in m["horizontes"] and "nublado" in m["horizontes"]["1"]
print("modelo.json:", os.path.getsize(os.path.join(tmp, "out", "modelo.json")), "B ·", "OK" if ok else "FALLO")
sys.exit(0 if ok else 1)
