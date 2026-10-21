"""Prueba de registro.py con cuadros reales ya publicados (carpeta con latest.json y f/).
Uso: python tests/test_registro.py DIR
Simula pasadas a las 11:20, 12:20, 13:20 y 14:20 y comprueba que se rellenan las etiquetas."""
import copy, gzip, json, os, sys, tempfile, time
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import registro as R, movimiento as M

d = sys.argv[1]
meta = json.load(open(os.path.join(d, "latest.json")))
G = meta["rejilla"]
cuadros = R.cargar_cuadros(meta, lambda f: open(os.path.join(d, f), "rb").read())
tmp = tempfile.mkdtemp()
for hhmm in ("11:20", "12:20", "13:20", "14:20"):
    corte = [c for c in meta["cuadros"] if c["t"][11:16] <= hhmm]
    m = copy.deepcopy(meta); m["cuadros"] = corte
    m["movimiento"] = M.campo_movimiento(corte, d, G["filas"], G["cols"], G, log=lambda *a: None)
    t = time.time(); r = R.registrar(m, {c["t"]: cuadros[c["t"]] for c in corte}, tmp)
    print(hhmm, "→", {k: r[k] for k in ("t0", "nueva", "etiquetas", "bytes", "puntos")}, f"{time.time()-t:.1f} s")
dia = [f for f in os.listdir(os.path.join(tmp, "registro")) if f.endswith(".gz")][0]
filas = R.leer_dia(os.path.join(tmp, "registro", dia))
print("filas del día:", len(filas), "| tamaño del archivo:", os.path.getsize(os.path.join(tmp, "registro", dia)), "B")
for f in filas:
    n = lambda k: sum(x is not None for x in f[k]) if f[k] else None
    print(f["t"], "f0:", n("f0"), "a1:", n("a1"), "a3:", n("a3"), "| etiquetas o1/o2/o3:", n("o1"), n("o2"), n("o3"))
# prueba de rotación
os.makedirs(os.path.join(tmp, "registro"), exist_ok=True)
for i in range(5):
    open(os.path.join(tmp, "registro", f"2026-01-0{i+1}.jsonl.gz"), "wb").write(os.urandom(1000))
b, tam = R.rotar(os.path.join(tmp, "registro"), max_bytes=3500)
print("rotación: borrados", b, "tamaño final", tam)
ok = filas[0]["o1"] is not None and len(b) >= 3
print("OK" if ok else "FALLO"); sys.exit(0 if ok else 1)
