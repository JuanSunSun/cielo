#!/usr/bin/env python3
"""Servidor estático de /data/_site + API mínima para el sitio de las alertas.

    GET  /api/alerta   -> {"soportado": true, "lat", "lon", "nombre", "radio"}  (lat null si no hay sitio)
    POST /api/alerta   -> guarda {"lat", "lon", "nombre", "radio"} en /data/alertas.json

Solo se acepta JSON (un formulario de otra web no puede cambiarlo) y se valida todo.
Tras guardar se crea /data/.reevaluar para que el bucle de alertas reaccione al momento.
"""
import functools
import http.server
import json
import os

RAIZ = os.environ.get("RAIZ", "/data/_site")
DATOS = os.path.dirname(RAIZ.rstrip("/")) or "/data"
RUTA_ALERTA = os.path.join(DATOS, "alertas.json")
PUERTO = int(os.environ.get("PUERTO", "8080"))


def leer_alerta():
    try:
        with open(RUTA_ALERTA, encoding="utf-8") as f:
            d = json.load(f)
    except (OSError, ValueError):
        d = {}
    return {"soportado": True, "lat": d.get("lat"), "lon": d.get("lon"),
            "nombre": d.get("nombre", ""), "radio": d.get("radio", 5)}


def validar(d):
    lat, lon = float(d["lat"]), float(d["lon"])
    if not (-90 <= lat <= 90 and -180 <= lon <= 180):
        raise ValueError("coordenadas fuera de rango")
    radio = float(d.get("radio", 5))
    if not (1 <= radio <= 50):
        raise ValueError("radio fuera de rango")
    nombre = str(d.get("nombre", "")).strip()[:60]
    return {"lat": round(lat, 4), "lon": round(lon, 4), "nombre": nombre, "radio": radio}


class Handler(http.server.SimpleHTTPRequestHandler):
    extensions_map = {
        **http.server.SimpleHTTPRequestHandler.extensions_map,
        ".webmanifest": "application/manifest+json",
        ".gz": "application/octet-stream",
    }

    def _json(self, codigo, obj):
        cuerpo = json.dumps(obj).encode()
        self.send_response(codigo)
        self.send_header("Content-Type", "application/json")
        self.send_header("Content-Length", str(len(cuerpo)))
        self.send_header("Cache-Control", "no-store")
        self.end_headers()
        self.wfile.write(cuerpo)

    def do_GET(self):
        if self.path.split("?")[0] == "/api/alerta":
            return self._json(200, leer_alerta())
        return super().do_GET()

    def do_POST(self):
        if self.path.split("?")[0] != "/api/alerta":
            return self._json(404, {"error": "no existe"})
        if "application/json" not in self.headers.get("Content-Type", ""):
            return self._json(415, {"error": "se espera JSON"})
        try:
            largo = int(self.headers.get("Content-Length", "0"))
            if not 0 < largo <= 2048:
                raise ValueError("tamaño")
            nuevo = validar(json.loads(self.rfile.read(largo)))
        except (ValueError, KeyError, TypeError):
            return self._json(400, {"error": "datos no válidos"})
        tmp = RUTA_ALERTA + ".tmp"
        with open(tmp, "w", encoding="utf-8") as f:
            json.dump(nuevo, f)
        os.replace(tmp, RUTA_ALERTA)
        open(os.path.join(DATOS, ".reevaluar"), "w").close()
        return self._json(200, dict(leer_alerta()))

    def end_headers(self):
        ruta = self.path.split("?")[0]
        if ruta.endswith("latest.json"):
            self.send_header("Cache-Control", "no-store")
        elif ruta.endswith("sw.js"):
            self.send_header("Cache-Control", "no-cache")
        super().end_headers()

    def log_message(self, *args):
        pass


if __name__ == "__main__":
    os.makedirs(RAIZ, exist_ok=True)
    handler = functools.partial(Handler, directory=RAIZ)
    http.server.ThreadingHTTPServer(("0.0.0.0", PUERTO), handler).serve_forever()
