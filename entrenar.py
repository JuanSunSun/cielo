"""Corrección estadística del pronóstico por advección.

Lee el registro (`registro/*.jsonl.gz`, generado por registro.py) y, para cada horizonte (+1, +2, +3 h),
ajusta una regresión logística que estima la probabilidad de que el punto esté
  · nublado   (> 50 % de nube en ~10 km)  y
  · despejado (≤ 10 %)
a partir de lo que sabía el pronóstico en ese momento. Compara con tres referencias (climatología,
persistencia y advección sin corregir) en un tramo final de días que el modelo no ha visto
(partición temporal, nunca aleatoria) y con un boosting de árboles como techo de referencia.

Solo escribe `modelo.json` (lo que carga la app) si hay datos suficientes y el modelo mejora de verdad
la advección (ya recalibrada con una sola variable) en el tramo de prueba (Brier). Si no, deja constancia en el informe y no toca nada.

Uso: python entrenar.py --registro datos/registro --salida datos
"""
import argparse
import gzip
import json
import math
import os
import sys
from datetime import datetime, timezone

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from registro import elev_solar, t_de  # noqa: E402

HORIZONTES = (1, 2, 3)
MIN_FILAS = 20000      # muestras con etiqueta por horizonte para entrenar
MIN_DIAS = 5           # días distintos en el registro
FRAC_PRUEBA = 0.25     # último tramo temporal reservado para evaluar
MEJORA_MIN = 0.03      # mejora relativa mínima del Brier frente a la advección sin corregir
OBJETIVOS = {"nublado": lambda o: o > 0.5, "despejado": lambda o: o <= 0.10}
NOMBRES = ["f0", "a", "cambio", "vel", "vel_x_cambio", "disp", "coh", "cerca", "sol", "lat", "lon_abs"]


def leer_registro(dir_reg):
    puntos = json.load(open(os.path.join(dir_reg, "puntos.json")))
    filas = []
    for f in sorted(os.listdir(dir_reg)):
        if f.endswith(".jsonl.gz"):
            with gzip.open(os.path.join(dir_reg, f), "rt", encoding="utf-8") as fh:
                filas += [json.loads(l) for l in fh if l.strip()]
    return puntos, filas


def construir(puntos, filas, h):
    """Matriz de características X, etiquetas observadas o (fracción), tiempos y días para el horizonte h."""
    lat = np.array([p[0] for p in puntos]); lon = np.array([p[1] for p in puntos])
    X, O, T = [], [], []
    for fila in filas:
        if fila.get(f"o{h}") is None:
            continue
        f0 = np.array([np.nan if x is None else x for x in fila["f0"]], float)
        a = np.array([np.nan if x is None else x for x in fila[f"a{h}"]], float)
        o = np.array([np.nan if x is None else x for x in fila[f"o{h}"]], float)
        u = np.array([np.nan if x is None else x for x in fila["u"]], float)
        v = np.array([np.nan if x is None else x for x in fila["v"]], float)
        d = np.array([20.0 if x is None else x for x in fila["d"]], float)
        k = np.array([0.0 if x is None else x for x in fila["k"]], float)
        c = np.array([200.0 if x is None else x for x in fila["c"]], float)
        ok = ~(np.isnan(f0) | np.isnan(a) | np.isnan(o) | np.isnan(u) | np.isnan(v))
        if not ok.any():
            continue
        t = t_de(fila["t"]); dt = (t.replace(tzinfo=timezone.utc))
        from datetime import timedelta
        sol = np.array([elev_solar(la, lo, dt + timedelta(hours=h)) for la, lo in zip(lat, lon)])
        vel = np.hypot(u, v)
        cambio = a - f0
        cols = np.column_stack([f0, a, cambio, vel, vel * cambio, d, k, np.minimum(c, 200.0), sol, lat, np.abs(lon)])
        X.append(cols[ok]); O.append(o[ok]); T += [fila["t"]] * int(ok.sum())
    if not X:
        return np.empty((0, len(NOMBRES))), np.empty(0), np.array([])
    return np.vstack(X), np.concatenate(O), np.array(T)


def brier(p, y):
    return float(np.mean((p - y) ** 2))


def logloss(p, y):
    p = np.clip(p, 1e-4, 1 - 1e-4)
    return float(-np.mean(y * np.log(p) + (1 - y) * np.log(1 - p)))


def entrenar_horizonte(X, O, T, h, informe):
    from sklearn.ensemble import HistGradientBoostingClassifier
    from sklearn.linear_model import LogisticRegression

    dias = sorted({t[:10] for t in T})
    corte_t = sorted(set(T))[int(len(set(T)) * (1 - FRAC_PRUEBA))]
    ent, pru = T < corte_t, T >= corte_t
    if ent.sum() < 100 or pru.sum() < 100:
        informe.append(f"+{h} h: pocos datos para partir en entrenamiento y prueba")
        return None
    media = X[ent].mean(0); desv = X[ent].std(0) + 1e-9
    Z = (X - media) / desv
    res = {"dias": len(dias), "n_entrenamiento": int(ent.sum()), "n_prueba": int(pru.sum()), "objetivos": {}}
    modelos = {}
    for nombre, f in OBJETIVOS.items():
        y = f(O).astype(int)
        if y[ent].min() == y[ent].max():
            continue
        clim = np.full(pru.sum(), y[ent].mean())
        # referencias honestas: persistencia y advección SIN corregir, pero calibradas con una logística
        # de una sola variable (así la mejora medida es la de añadir información, no la de recalibrar)
        def una_variable(col):
            m1 = LogisticRegression(C=10.0, max_iter=300).fit(X[ent][:, [col]], y[ent])
            return m1.predict_proba(X[pru][:, [col]])[:, 1]
        pers_p, adv_p = una_variable(0), una_variable(1)
        lr = LogisticRegression(C=1.0, max_iter=500).fit(Z[ent], y[ent])
        p_lr = lr.predict_proba(Z[pru])[:, 1]
        gb = HistGradientBoostingClassifier(max_iter=150, learning_rate=0.08, max_leaf_nodes=15).fit(X[ent], y[ent])
        p_gb = gb.predict_proba(X[pru])[:, 1]
        m = {n: {"brier": brier(p, y[pru]), "logloss": logloss(p, y[pru])} for n, p in
             (("climatología", clim), ("persistencia", pers_p), ("advección", adv_p),
              ("logística", p_lr), ("boosting", p_gb))}
        m["mejora_vs_advección"] = 1 - m["logística"]["brier"] / m["advección"]["brier"]
        res["objetivos"][nombre] = m
        modelos[nombre] = {"intercepto": float(lr.intercept_[0]), "coef": [float(c) for c in lr.coef_[0]],
                           "mejora": float(m["mejora_vs_advección"])}
    res["media"] = [float(x) for x in media]; res["desv"] = [float(x) for x in desv]
    informe.append(f"+{h} h  ({res['n_entrenamiento']} entrenamiento / {res['n_prueba']} prueba, {len(dias)} días)")
    for nombre, m in res["objetivos"].items():
        fila = "   ".join(f"{k} {m[k]['brier']:.4f}" for k in ("climatología", "persistencia", "advección", "logística", "boosting"))
        informe.append(f"   {nombre:9s} Brier: {fila}   → mejora logística vs advección {m['mejora_vs_advección']*100:+.1f} %")
    return res, modelos


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--registro", default="datos/registro")
    ap.add_argument("--salida", default="datos")
    ap.add_argument("--forzar", action="store_true", help="ignora los mínimos de datos (solo para pruebas)")
    a = ap.parse_args(argv)
    puntos, filas = leer_registro(a.registro)
    dias = sorted({f["t"][:10] for f in filas})
    informe = [f"Registro: {len(filas)} pasadas, {len(dias)} días ({dias[0] if dias else '-'} a {dias[-1] if dias else '-'})"]
    modelo = {"version": 1, "entrenado": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%MZ"),
              "caracteristicas": NOMBRES, "horizontes": {}}
    publicable = False
    for h in HORIZONTES:
        X, O, T = construir(puntos, filas, h)
        informe.append("")
        if (len(O) < MIN_FILAS or len(dias) < MIN_DIAS) and not a.forzar:
            informe.append(f"+{h} h: {len(O)} muestras con etiqueta y {len(dias)} días; hacen falta {MIN_FILAS} y {MIN_DIAS}. Sin entrenar.")
            continue
        r = entrenar_horizonte(X, O, T, h, informe)
        if not r:
            continue
        res, modelos = r
        buenos = {n: m for n, m in modelos.items() if m["mejora"] >= MEJORA_MIN}
        if buenos:
            modelo["horizontes"][str(h)] = {"media": res["media"], "desv": res["desv"], **buenos}
            publicable = True
        for n, m in modelos.items():
            if m["mejora"] < MEJORA_MIN:
                informe.append(f"   {n}: mejora {m['mejora']*100:+.1f} % < {MEJORA_MIN*100:.0f} %, no se publica")
    texto = "\n".join(informe)
    print(texto)
    os.makedirs(a.salida, exist_ok=True)
    open(os.path.join(a.salida, "informe.txt"), "w", encoding="utf-8").write(texto + "\n")
    ruta = os.path.join(a.salida, "modelo.json")
    if publicable:
        json.dump(modelo, open(ruta, "w"), separators=(",", ":"))
        print("→ modelo.json escrito")
    elif os.path.exists(ruta):
        print("→ se conserva el modelo.json anterior")
    return 0


if __name__ == "__main__":
    sys.exit(main())
