#!/bin/sh
# Un contenedor: un servidor Python sirve /data/_site y un bucle regenera los datos cada 15 min.
set -u
SALIDA=/data/_site
INTERVALO="${INTERVALO_SEG:-900}"
mkdir -p "$SALIDA"

publicar_app() { cp -r /opt/cielo/app/. "$SALIDA"/; }

ciclo() {
  sleep 2
  while true; do
    echo "[cielo] $(date -u +%FT%TZ) procesando GOES-19"
    python /opt/cielo/procesar.py --salida "$SALIDA" --sitio "http://127.0.0.1:8080" \
      || echo "[cielo] procesar.py falló; reintento en $INTERVALO s"
    publicar_app
    python /opt/cielo/alertas.py --salida "$SALIDA" || echo "[cielo] alertas.py falló"
    esperar
  done
}

# espera INTERVALO s, pero si la app cambia el sitio de las alertas (/data/.reevaluar) las recalcula al momento
esperar() {
  n=0
  while [ "$n" -lt "$INTERVALO" ]; do
    if [ -f /data/.reevaluar ]; then
      rm -f /data/.reevaluar
      python /opt/cielo/alertas.py --salida "$SALIDA" || echo "[cielo] alertas.py falló"
    fi
    sleep 5; n=$((n + 5))
  done
}

publicar_app          # la web carga aunque aún no haya datos
ciclo &
exec python /opt/cielo/servidor.py
