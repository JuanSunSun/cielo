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

## Movimiento de las nubes y pronóstico a corto plazo

`movimiento.py` mide hacia dónde se mueven las nubes comparando el cuadro más reciente con el de ~30 min antes
(correlación de fase por teselas de ~210 km; los bordes rectos, que solo fijan la componente perpendicular,
se resuelven con las teselas vecinas por mínimos cuadrados robustos). El resultado va en `data/latest.json`
(`movimiento.vectores`, en km/h este/norte).

La app, para el punto elegido, mira el área aguas arriba (`posición − vector·t`, con el radio ensanchado un 20 %
de la distancia recorrida) cada 15 min hasta 3,5 h y resume cuándo entran o se van las nubes. Es una extrapolación:
no prevé nubes que se formen o se disipen, ni las ligadas al relieve. Fiable hasta ~2 h.

Prueba del módulo: `python tests/test_movimiento.py`.

### Cuánto se puede fiar (pruebas retrospectivas con GOES-19 real)

`python tests/hindcast.py DIR` repite el pronóstico sobre cuadros ya publicados y lo compara con lo que ocurrió
(DIR = carpeta con `latest.json` y `f/`; el flujo «Diagnóstico» los copia a la rama `diag-out`).
Primera medición, 3 h de datos del 9-oct-2026 a mediodía, Chile completo, fracción de nube en ~20 km:

| Horizonte | Error de la fracción: «no cambia» → advección | Punto despejado que pasa a nublado: detectado / falsas alarmas |
|---|---|---|
| +30 min | 7,3 % → 6,7 % | 25 % / 67 % |
| +60 min | 11,2 % → 10,5 % | 27 % / 63 % |
| +90 min | 13,9 % → 13,2 % | 25 % / 66 % |

Con nubes que se mueven > 10 km/h y movimiento repetido entre parejas de imágenes sube a ~60 % detectado
(con ~65 % de falsas alarmas, porque el aviso salta a partir del 25 % previsto). La mayor parte de la región
(costa, cordillera) tiene nubes casi estacionarias o que cambian de forma, donde extrapolar apenas mejora.
La app marca «Orientativo» cuando no se da esa condición.

Pendiente para mejorar: seguir el movimiento con el infrarrojo (BT 11 µm) en vez de la máscara binaria,
acumular semanas de verificación para calibrar probabilidades y mezclar con modelos numéricos más allá de 2 h.
