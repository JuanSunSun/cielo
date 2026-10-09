# Cielo

¿Está despejado? Web app instalable (Chrome en Mac, iPhone) que muestra las nubes de
GOES-19 sobre Chile y la Patagonia **en su posición real** (corregidas de paralaje) y da
un veredicto para cualquier punto del mapa.

**App:** https://juansunsun.github.io/cielo/

## Cómo funciona

```
GitHub Actions (cada 15 min)                    App (navegador)
  procesar.py                                     index.html
  ├─ GOES-19 ABI-L2-ACMF  máscara de nubes        ├─ descarga data/latest.json + cuadros
  ├─ ABI-L2-CMIPF C14/C07 11.2 y 3.9 µm           ├─ pinta la capa de nubes (Leaflet)
  ├─ ABI-L2-ACHAF         altura del tope         ├─ clic → % de nubes en el radio,
  ├─ paralaje: nube → posición real               │   nube más cercana, velocidad y ETA
  ├─ niebla nocturna: BT11 − BT3.9 > 2.5 K        ├─ Open-Meteo: nubes por capas, rocío
  └─ rejilla 0.02° → data/f/*.bin.gz  ──Pages──▶  └─ modo rojo, sitios guardados, offline
```

Códigos de la rejilla: 0 sin dato · 1 despejado · 2 despejado dudoso (prob. de nube 0.35–0.5) · 3 nube baja
(<2.5 km) · 4 media · 5 alta (>6 km) · 6 niebla/estrato nocturno.

El histórico (~4 h) no se guarda en git: cada ejecución lo recupera del sitio publicado.

## Instalar la app

- **Mac (Chrome):** abre la URL → icono de instalar en la barra de direcciones (o menú ⋮ → *Transmitir, guardar y compartir* → *Instalar página como app*).
- **iPhone (Safari o Chrome):** abre la URL → Compartir → *Añadir a pantalla de inicio*.

## Límites

- Cobertura 17.5°S–56°S, 64°W–76°W. Euskadi no: GOES-19 la ve casi rasante; haría falta MTG/Meteosat (EUMETSAT, con claves).
- Latencia: imagen de ~10–20 min + hasta 15 min de la tarea programada (GitHub puede retrasarla).
- Resolución ~2 km. El test de niebla puede dar falsos positivos sobre desiertos (Atacama) por la emisividad de la arena a 3.9 µm.
- GitHub desactiva las tareas programadas tras 60 días sin commits: si pasa, *Actions → Mapa de nubes → Enable workflow*.

## Desarrollo local

```
pip install -r requirements.txt
python procesar.py --salida _site      # descarga GOES-19 real
cp -r app/. _site/ && python -m http.server -d _site 8000
```
