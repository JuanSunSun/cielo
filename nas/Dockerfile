FROM python:3.12-slim

WORKDIR /opt/cielo
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

COPY procesar.py .
COPY movimiento.py .
COPY alertas.py .
COPY app ./app
COPY nas/servidor.py .
COPY nas/entrypoint.sh /entrypoint.sh
RUN chmod +x /entrypoint.sh

ENV TZ=America/Santiago
EXPOSE 8080
ENTRYPOINT ["/entrypoint.sh"]
