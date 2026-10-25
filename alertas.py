#!/usr/bin/env python3
"""
alertas.py — Avisos de Cielo por Telegram para UN sitio.

Se ejecuta tras cada ciclo de procesar.py (cada 15 min). Lee la última rejilla de
nubes, estima con el viento de cada altura (Open-Meteo, 850/700/300 hPa) qué nubes
llegarán al sitio en las próximas horas —la misma lógica que el panel de la app— y
envía un mensaje solo cuando cambia algo:

    despejar   ahora hay nubes y en 1–3 h se despeja
    nublar     ahora está despejado y en 1–3 h llegan nubes
    niebla     el satélite ve niebla/estrato nocturno sobre el sitio
    rocio      temperatura y punto de rocío casi iguales (riesgo en la óptica)
    lluvia     ☔ lluvia prevista en las próximas horas (o ya cayendo); a cualquier hora
    niebla_local  la estación local marca condiciones de niebla y el satélite no la ve

Si /estacion/latest.json (receptor Ecowitt en el NAS) tiene menos de 10 min, rocío y lluvia usan las
medidas de la estación en vez del modelo, y cada ciclo añade una fila a /data/comparacion_estacion.csv
(satélite vs estación) para recalibrar umbrales. Instalación en tejado, no estándar: valores válidos
para tendencias y umbrales locales, no homologables.

Cada tipo tiene un enfriamiento para no repetirse. Por defecto despejar/nublar/niebla/rocío
solo avisan de noche: el evento debe ocurrir con el sol a menos de -6°, y el aviso se abre
2 h antes del atardecer (mirando hasta 6 h hacia delante) para poder organizar la tarde.
La lluvia no depende de la hora.

El sitio se cambia desde la app (botón "Avisarme de este sitio"), que lo guarda en
/data/alertas.json. Si ese archivo no existe se usan ALERTA_* del entorno.

Configuración por variables de entorno (se ponen en alertas.env, fuera de git):
    ALERTA_LAT, ALERTA_LON        sitio por defecto (si no hay alertas.json)
    ALERTA_NOMBRE                 nombre en los mensajes (por defecto "tu sitio")
    ALERTA_RADIO                  radio de análisis en km (por defecto 5)
    TELEGRAM_TOKEN, TELEGRAM_CHAT_ID
    CIELO_URL                     enlace a la app, opcional
    SOLO_NOCHE                    1 (por defecto) o 0
    ANTES_ATARDECER_H             horas antes del atardecer en que se abre el aviso (por defecto 2)
    LLUVIA_HORAS                  horizonte del aviso de lluvia (por defecto 3)

Uso:
    python alertas.py --salida /data/_site          ciclo normal
    python alertas.py --salida DIR --prueba         muestra qué enviaría, sin enviar ni guardar
    python alertas.py --chat-id                     imprime tu chat_id (tras escribirle al bot)
    python alertas.py --enviar-prueba               manda "avisos activados" para comprobar
"""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import os
import sys
from datetime import datetime, timedelta, timezone

import numpy as np
import requests

NIVELES = [  # (hPa, clases de la rejilla, nombre)
    (850, (3, 6), "bajas"),
    (700, (4,), "medias"),
    (300, (5,), "altas"),
]
NUBE = (3, 4, 5, 6)
COOLDOWN_H = {"despejar": 4, "nublar": 3, "niebla": 6, "niebla_local": 6, "rocio": 6, "lluvia": 6}
ESTACION_JSON = os.environ.get("ESTACION_JSON", "/estacion/latest.json")
ESTACION_MAX_EDAD_S = 600
VENTANA_H = 3.5          # horas hacia delante que se consideran para avisar
VENTANA_TARDE_H = 6.0    # ... antes de que anochezca, para planificar la noche


def log(*a):
    print("[alertas]", *a, flush=True)


# ───────────────────────── datos ─────────────────────────
def cargar_ultimo(salida):
    with open(os.path.join(salida, "data", "latest.json"), encoding="utf-8") as f:
        meta = json.load(f)
    if not meta.get("cuadros"):
        return None
    ult = meta["cuadros"][-1]
    with open(os.path.join(salida, "data", ult["f"]), "rb") as f:
        raw = f.read()
    if raw[:2] == b"\x1f\x8b":
        raw = gzip.decompress(raw)
    G = meta["rejilla"]
    g = np.frombuffer(raw, dtype=np.uint8).reshape(G["filas"], G["cols"])
    t = datetime.strptime(ult["t"], "%Y-%m-%dT%H:%MZ").replace(tzinfo=timezone.utc)
    return g, G, t


def pedir_modelo(lat, lon):
    url = (f"https://api.open-meteo.com/v1/forecast?latitude={lat:.3f}&longitude={lon:.3f}"
           "&current=temperature_2m,dew_point_2m,precipitation"
           "&hourly=wind_speed_850hPa,wind_direction_850hPa,wind_speed_700hPa,wind_direction_700hPa,"
           "wind_speed_300hPa,wind_direction_300hPa,precipitation,precipitation_probability"
           "&forecast_hours=12&timezone=auto")
    r = requests.get(url, timeout=20)
    r.raise_for_status()
    return r.json()


def leer_estacion(ahora):
    """Última lectura de la estación local (dict de la vista obs), o None si falta o tiene >10 min."""
    try:
        with open(ESTACION_JSON, encoding="utf-8") as f:
            o = json.load(f)["obs"]
        t = datetime.strptime(o["time_utc"], "%Y-%m-%dT%H:%M:%SZ").replace(tzinfo=timezone.utc)
    except OSError:
        return None
    except (ValueError, KeyError, TypeError) as e:
        log(f"latest.json de la estación ilegible ({e}): se usa solo el modelo.")
        return None
    edad = (ahora - t).total_seconds()
    if edad < -60 or edad > ESTACION_MAX_EDAD_S:
        log(f"Estación sin datos recientes (última hace {int(edad)} s): se usa solo el modelo.")
        return None
    o["edad_s"] = int(edad)
    return o


def resumen_estacion(e):
    partes = []
    if e.get("temp_c") is not None:
        partes.append(f"{e['temp_c']:.1f} °C")
    if e.get("humidity_pct") is not None:
        partes.append(f"HR {e['humidity_pct']:.0f} %")
    if e.get("dewpoint_c") is not None:
        partes.append(f"rocío {e['dewpoint_c']:.1f} °C")
    if e.get("wind_ms") is not None:
        partes.append(f"viento {e['wind_ms']:.1f} m/s")
    return "Estación: " + ", ".join(partes) + "."


PARES_COLS = ["ts_utc", "frame_utc", "sat_nube_frac", "sat_niebla_frac", "sat_despejado_frac",
              "st_time_utc", "st_edad_s", "st_temp_c", "st_humidity_pct", "st_dewpoint_c", "st_spread_c",
              "st_wind_ms", "st_gust_ms", "st_solar_wm2", "st_rain_rate_mm_h", "st_fog_flag",
              "lat", "lon", "radio_km"]


def registrar_pares(ruta, ahora, t_frame, g, G, cfg, estacion):
    """Añade (satélite, estación) al CSV para calibrar umbrales más adelante.
    Las fracciones del satélite son del radio configurado; vacío = fuera de cobertura o sin dato."""
    def sat(clases):
        v = muestrear(g, G, cfg["lat"], cfg["lon"], cfg["radio"], clases)
        return "" if v is None else round(v, 3)

    e = estacion or {}

    def st(k):
        v = e.get(k)
        return "" if v is None else v

    fila = [ahora.strftime("%Y-%m-%dT%H:%M:%SZ"), t_frame.strftime("%Y-%m-%dT%H:%MZ"),
            sat(NUBE), sat((6,)), sat((1, 2)),
            st("time_utc"), st("edad_s"), st("temp_c"), st("humidity_pct"), st("dewpoint_c"), st("spread_c"),
            st("wind_ms"), st("gust_ms"), st("solar_wm2"), st("rain_rate_mm_h"), st("fog_flag"),
            cfg["lat"], cfg["lon"], cfg["radio"]]
    nuevo = not os.path.exists(ruta)
    with open(ruta, "a", newline="", encoding="utf-8") as f:
        if nuevo:
            f.write("# Estación en tejado, no estándar; valores válidos para tendencias y umbrales locales, "
                    "no homologables. Satélite = fracción de celdas en el radio (clases 3-6 nube, 6 niebla, 1-2 despejado).\n")
        w = csv.writer(f)
        if nuevo:
            w.writerow(PARES_COLS)
        w.writerow(fila)


# ───────────────────────── cálculo ─────────────────────────
def elev_solar(lat, lon, t):
    """Misma fórmula que elevSolar() de la app. t: datetime con tz UTC."""
    doy = t.timetuple().tm_yday
    h = t.hour + t.minute / 60
    g = 2 * math.pi / 365 * (doy - 1 + (h - 12) / 24)
    dec = (.006918 - .399912 * math.cos(g) + .070257 * math.sin(g) - .006758 * math.cos(2 * g)
           + .000907 * math.sin(2 * g) - .002697 * math.cos(3 * g) + .00148 * math.sin(3 * g))
    eq = 229.18 * (.000075 + .001868 * math.cos(g) - .032077 * math.sin(g)
                   - .014615 * math.cos(2 * g) - .040849 * math.sin(2 * g))
    ha = math.radians(((h * 60 + eq + 4 * lon) / 4) - 180)
    la = math.radians(lat)
    return 90 - math.degrees(math.acos(math.sin(la) * math.sin(dec) + math.cos(la) * math.cos(dec) * math.cos(ha)))


def ventana_aviso(lat, lon, ahora, horas_antes):
    """True si ya anocheció (sol bajo el horizonte) o faltan <= horas_antes para la puesta de sol."""
    if elev_solar(lat, lon, ahora) < 0:
        return True
    for k in range(1, 24 * 12 + 1):                      # busca la puesta de sol en pasos de 5 min
        t = ahora + timedelta(minutes=5 * k)
        if elev_solar(lat, lon, t) < 0:
            return (t - ahora).total_seconds() / 3600 <= horas_antes
    return False                                          # el sol no se pone (verano polar)


def muestrear(g, G, lat, lon, radio, clases):
    """Fracción de celdas válidas de `clases` en un círculo. None si el centro cae fuera."""
    paso = G["paso"]
    kx = 111.32 * math.cos(math.radians(lat)) * paso
    ky = 110.57 * paso
    r0 = (G["lat_max"] - lat) / paso - .5
    c0 = (lon - G["lon_min"]) / paso - .5
    if r0 < 0 or r0 >= G["filas"] or c0 < 0 or c0 >= G["cols"]:
        return None
    dr, dc = math.ceil(radio / ky), math.ceil(radio / kx)
    r1, r2 = max(0, math.floor(r0 - dr)), min(G["filas"] - 1, math.ceil(r0 + dr))
    c1, c2 = max(0, math.floor(c0 - dc)), min(G["cols"] - 1, math.ceil(c0 + dc))
    sub = g[r1:r2 + 1, c1:c2 + 1]
    rr, cc = np.mgrid[r1:r2 + 1, c1:c2 + 1]
    dentro = np.hypot((cc - c0) * kx, (r0 - rr) * ky) <= radio
    valido = (sub > 0) & dentro
    n = int(valido.sum())
    if n == 0:  # radio menor que una celda: celda más cercana
        v = int(g[min(G["filas"] - 1, max(0, round(r0))), min(G["cols"] - 1, max(0, round(c0)))])
        return None if v == 0 else (1.0 if v in clases else 0.0)
    return float(np.isin(sub[valido], clases).sum()) / n


def viento_hacia(h, nivel, n):
    """Vector (este, norte) km/h hacia DE DONDE viene el viento, media de las horas 0..n."""
    sp_l, d_l = h.get(f"wind_speed_{nivel}hPa"), h.get(f"wind_direction_{nivel}hPa")
    if not sp_l or not d_l:
        return None
    e = no = k = 0
    for j in range(min(n, len(sp_l) - 1) + 1):
        sp, d = sp_l[j], d_l[j]
        if sp is None or d is None:
            continue
        e += sp * math.sin(math.radians(d))
        no += sp * math.cos(math.radians(d))
        k += 1
    return (e / k, no / k) if k else None


def pronostico(g, G, t_frame, modelo, lat, lon, radio):
    """Lista de {hora_local, t_utc, frac} para las próximas horas (frac None = sin datos)."""
    h = modelo["hourly"]
    off = modelo.get("utc_offset_seconds", 0)
    filas = []
    for i, tl in enumerate(h["time"]):
        if len(filas) >= 6:
            break
        t_utc = datetime.strptime(tl, "%Y-%m-%dT%H:%M").replace(tzinfo=timezone.utc) - timedelta(seconds=off)
        dt = (t_utc - t_frame).total_seconds() / 3600
        if dt < 0.4:
            continue
        libres = []
        for nivel, clases, _ in NIVELES:
            w = viento_hacia(h, nivel, i)
            if not w:
                continue
            lat2 = lat + w[1] * dt / 110.57
            lon2 = lon + w[0] * dt / (111.32 * math.cos(math.radians(lat)))
            f = muestrear(g, G, lat2, lon2, radio, clases)
            if f is not None:
                libres.append(1 - f)
        frac = None if not libres else 1 - math.prod(libres)
        filas.append({"hora": tl[11:16], "t": t_utc, "frac": frac})
    return filas


# ───────────────────────── decisión ─────────────────────────
def decidir(ahora, t_frame, g, G, modelo, cfg, estado, estacion=None):
    lat, lon, radio, nombre = cfg["lat"], cfg["lon"], cfg["radio"], cfg["nombre"]
    a_nube = muestrear(g, G, lat, lon, radio, NUBE)
    if a_nube is None:
        log("El sitio está fuera de la cobertura.")
        return []
    a_niebla = muestrear(g, G, lat, lon, radio, (6,)) or 0.0
    antes_noche = elev_solar(lat, lon, ahora) > -6
    horizonte = VENTANA_TARDE_H if (cfg["solo_noche"] and antes_noche) else VENTANA_H
    filas = [f for f in pronostico(g, G, t_frame, modelo, lat, lon, radio)
             if f["frac"] is not None and (f["t"] - ahora).total_seconds() / 3600 <= horizonte]
    en_ventana = (not cfg["solo_noche"]) or ventana_aviso(lat, lon, ahora, cfg["antes_atardecer_h"])
    msgs = []

    def noche_en(t):                       # niebla y rocío: es de noche en ese instante
        return (not cfg["solo_noche"]) or elev_solar(lat, lon, t) < -6

    def evento_ok(t):                      # despejar / nublar: ventana abierta y el evento ocurre de noche
        return en_ventana and noche_en(t)

    def libre(tipo):
        ult = estado.get("ultimo", {}).get(tipo)
        return not ult or (ahora - datetime.fromisoformat(ult)).total_seconds() / 3600 >= COOLDOWN_H[tipo]

    # se va a despejar: ahora nublado y el pronóstico llega a despejado y no vuelve a nublarse
    if a_nube > .5 and libre("despejar"):
        for i, f in enumerate(filas):
            resto = [x["frac"] for x in filas[i:]]
            if f["frac"] <= .10 and max(resto) <= .5 and evento_ok(f["t"]):
                msgs.append(("despejar", f"☀️ Se despeja sobre {nombre} hacia las {f['hora']}: "
                             f"~{round(f['frac'] * 100)} % de nubes (ahora {round(a_nube * 100)} %)."))
                break

    # se va a nublar: ahora despejado y llegan nubes
    if a_nube <= .25 and libre("nublar"):
        for f in filas:
            if f["frac"] > .5 and evento_ok(f["t"]):
                msgs.append(("nublar", f"☁️ Llegan nubes a {nombre} hacia las {f['hora']}: "
                             f"~{round(f['frac'] * 100)} % (ahora {round(a_nube * 100)} %)."))
                break

    # niebla / estrato nocturno visto por satélite
    if a_niebla >= .3 and libre("niebla") and noche_en(ahora):
        texto = (f"🌫️ El satélite ve niebla o estrato bajo sobre {nombre} "
                 f"({round(a_niebla * 100)} % del radio).")
        if estacion:
            texto += " " + resumen_estacion(estacion)
        msgs.append(("niebla", texto))
    # niebla en la estación que el satélite no ve (nubes altas encima, estrato fino...)
    elif estacion and estacion.get("fog_flag") == 1 and a_niebla < .3 and libre("niebla_local") \
            and noche_en(ahora):
        msgs.append(("niebla_local", f"🌫️ Condiciones de niebla en la estación de {nombre}. "
                     f"{resumen_estacion(estacion)} El satélite no la ve ({round(a_niebla * 100)} % del radio)."))

    # ☔ lluvia: ya cayendo, o prevista en las próximas horas (no depende de la hora del día)
    if libre("lluvia"):
        cur = modelo.get("current") or {}
        h = modelo["hourly"]
        off = modelo.get("utc_offset_seconds", 0)
        ahora_mm = cur.get("precipitation") or 0
        llueve_est = estacion is not None and (estacion.get("rain_rate_mm_h") or 0) >= 0.1
        if llueve_est:
            msgs.append(("lluvia", f"☔ LLUVIA en {nombre}: la estación mide {estacion['rain_rate_mm_h']:g} mm/h "
                         "ahora. Protege el equipo."))
        elif ahora_mm >= 0.1:
            msgs.append(("lluvia", f"☔ LLUVIA en {nombre}: está cayendo ahora ({ahora_mm:g} mm). Protege el equipo."))
        else:
            for tl, mm, pr in zip(h["time"], h.get("precipitation") or [], h.get("precipitation_probability") or []):
                t_utc = datetime.strptime(tl, "%Y-%m-%dT%H:%M").replace(tzinfo=timezone.utc) - timedelta(seconds=off)
                dt = (t_utc - ahora).total_seconds() / 3600
                if dt < -0.9 or dt > cfg["lluvia_horas"] or mm is None:
                    continue
                pr = pr if pr is not None else 100
                if (mm >= 0.3 and pr >= 40) or (mm >= 0.2 and pr >= 70):
                    msgs.append(("lluvia", f"☔ LLUVIA prevista en {nombre} hacia las {tl[11:16]}: "
                                 f"~{mm:g} mm/h, probabilidad {pr} %. Protege el equipo."))
                    break

    # riesgo de rocío
    cur = modelo.get("current") or {}
    t_aire, t_rocio, fuente = cur.get("temperature_2m"), cur.get("dew_point_2m"), ""
    if estacion and estacion.get("temp_c") is not None and estacion.get("dewpoint_c") is not None:
        t_aire, t_rocio, fuente = estacion["temp_c"], estacion["dewpoint_c"], " (estación)"
    if t_aire is not None and t_rocio is not None and t_aire - t_rocio < 2 and libre("rocio") \
            and noche_en(ahora):
        msgs.append(("rocio", f"💧 Riesgo de rocío en {nombre}: {t_aire:g} °C con punto de "
                     f"rocío {t_rocio:g} °C{fuente}. Ojo con la óptica."))
    return msgs


# ───────────────────────── Telegram ─────────────────────────
def enviar(token, chat_id, texto):
    r = requests.post(f"https://api.telegram.org/bot{token}/sendMessage",
                      json={"chat_id": chat_id, "text": texto, "disable_web_page_preview": True}, timeout=20)
    if not r.ok:
        # nunca imprimir la URL (lleva el token)
        raise RuntimeError(f"Telegram respondió {r.status_code}: {r.text[:200]}")


def chat_ids(token):
    r = requests.get(f"https://api.telegram.org/bot{token}/getUpdates", timeout=20)
    r.raise_for_status()
    vistos = {}
    for u in r.json().get("result", []):
        m = u.get("message") or u.get("edited_message") or {}
        c = m.get("chat")
        if c:
            vistos[c["id"]] = c.get("first_name") or c.get("title") or c.get("username") or "?"
    return vistos


# ───────────────────────── principal ─────────────────────────
def leer_config(ruta_json=None):
    """El sitio sale de alertas.json (lo cambia la app); si no existe, de ALERTA_* del entorno."""
    sitio = {}
    if ruta_json:
        try:
            with open(ruta_json, encoding="utf-8") as f:
                sitio = json.load(f)
        except (OSError, ValueError):
            sitio = {}
    try:
        lat = float(sitio.get("lat", os.environ.get("ALERTA_LAT")))
        lon = float(sitio.get("lon", os.environ.get("ALERTA_LON")))
    except (TypeError, ValueError):
        return None
    return {
        "lat": lat, "lon": lon,
        "nombre": sitio.get("nombre") or os.environ.get("ALERTA_NOMBRE", "tu sitio"),
        "radio": float(sitio.get("radio") or os.environ.get("ALERTA_RADIO", "5")),
        "token": os.environ.get("TELEGRAM_TOKEN", ""),
        "chat": os.environ.get("TELEGRAM_CHAT_ID", ""),
        "url": os.environ.get("CIELO_URL", ""),
        "solo_noche": os.environ.get("SOLO_NOCHE", "1") != "0",
        "antes_atardecer_h": float(os.environ.get("ANTES_ATARDECER_H", "2")),
        "lluvia_horas": float(os.environ.get("LLUVIA_HORAS", "3")),
    }


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--salida", default="/data/_site")
    ap.add_argument("--estado", default=None, help="archivo de estado (por defecto <salida>/../alertas_estado.json)")
    ap.add_argument("--prueba", action="store_true", help="no envía ni guarda: solo muestra")
    ap.add_argument("--chat-id", action="store_true")
    ap.add_argument("--enviar-prueba", action="store_true")
    args = ap.parse_args(argv)

    token = os.environ.get("TELEGRAM_TOKEN", "")
    if args.chat_id:
        if not token:
            sys.exit("Falta TELEGRAM_TOKEN en el entorno.")
        ids = chat_ids(token)
        if not ids:
            sys.exit("No hay mensajes: escríbele cualquier cosa a tu bot en Telegram y repite.")
        for cid, nom in ids.items():
            print(f"chat_id={cid}  ({nom})")
        return

    ruta_json = os.path.join(os.path.dirname(os.path.abspath(args.salida)), "alertas.json")
    cfg = leer_config(ruta_json)
    if not cfg:
        log("Sin sitio (alertas.json o ALERTA_LAT/ALERTA_LON): alertas desactivadas.")
        return
    if args.enviar_prueba:
        enviar(cfg["token"], cfg["chat"], f"✅ Cielo: avisos activados para {cfg['nombre']}.")
        log("Mensaje de prueba enviado.")
        return
    if not args.prueba and not (cfg["token"] and cfg["chat"]):
        log("Faltan TELEGRAM_TOKEN / TELEGRAM_CHAT_ID: alertas desactivadas.")
        return

    ruta_estado = args.estado or os.path.join(os.path.dirname(os.path.abspath(args.salida)), "alertas_estado.json")
    try:
        with open(ruta_estado, encoding="utf-8") as f:
            estado = json.load(f)
    except (OSError, ValueError):
        estado = {}

    sitio_id = f"{cfg['lat']:.4f},{cfg['lon']:.4f},{cfg['radio']:g}"
    if estado.get("sitio") != sitio_id:
        estado = {"sitio": sitio_id, "ultimo": {}}
        texto = f"📍 Avisos de Cielo ahora para {cfg['nombre']} ({cfg['lat']:.3f}, {cfg['lon']:.3f}), radio {cfg['radio']:g} km."
        if args.prueba:
            log(texto)
        else:
            try:
                enviar(cfg["token"], cfg["chat"], texto)
            except Exception as e:  # noqa: BLE001
                log(f"No se pudo enviar el aviso de cambio de sitio: {e}")

    try:
        ult = cargar_ultimo(args.salida)
    except (OSError, ValueError, KeyError) as e:
        log(f"Aún no hay datos ({e}).")
        return
    if not ult:
        return
    g, G, t_frame = ult
    ahora = datetime.now(timezone.utc)
    if (ahora - t_frame) > timedelta(hours=2):
        log(f"La última imagen es de hace {(ahora - t_frame)}; no se avisa con datos viejos.")
        return
    estacion = leer_estacion(ahora)
    if not args.prueba:
        ruta_pares = os.path.join(os.path.dirname(os.path.abspath(args.salida)), "comparacion_estacion.csv")
        try:
            registrar_pares(ruta_pares, ahora, t_frame, g, G, cfg, estacion)
        except Exception as e:  # noqa: BLE001
            log(f"No se pudo registrar el par satélite/estación: {e}")
    modelo = pedir_modelo(cfg["lat"], cfg["lon"])
    msgs = decidir(ahora, t_frame, g, G, modelo, cfg, estado, estacion)
    if not msgs:
        log("Nada que avisar.")
    for tipo, texto in msgs:
        if cfg["url"]:
            texto += f"\n{cfg['url']}"
        if args.prueba:
            log(f"[{tipo}] {texto}")
            continue
        try:
            enviar(cfg["token"], cfg["chat"], texto)
            estado.setdefault("ultimo", {})[tipo] = ahora.isoformat()
            log(f"Enviado: {tipo}")
        except Exception as e:  # noqa: BLE001
            log(f"No se pudo enviar ({tipo}): {e}")
    if not args.prueba:
        with open(ruta_estado, "w", encoding="utf-8") as f:
            json.dump(estado, f)


if __name__ == "__main__":
    main()
