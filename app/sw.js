// Service worker de Cielo: la app funciona sin red con los últimos datos vistos.
const VERSION = "cielo-v3";
const BASE = ["./", "index.html", "manifest.webmanifest", "icono-192.png", "icono-512.png", "icono-180.png"];

self.addEventListener("install", e => {
  e.waitUntil(caches.open(VERSION).then(c => c.addAll(BASE)).then(() => self.skipWaiting()));
});
self.addEventListener("activate", e => {
  e.waitUntil(caches.keys().then(ks => Promise.all(ks.filter(k => k !== VERSION).map(k => caches.delete(k))))
    .then(() => self.clients.claim()));
});
self.addEventListener("fetch", e => {
  const url = new URL(e.request.url);
  if (e.request.method !== "GET") return;
  // Cuadros de nubes: nombres únicos, nunca cambian -> primero caché
  if (url.pathname.includes("/data/f/")) {
    e.respondWith(caches.open(VERSION).then(async c => {
      const hit = await c.match(e.request);
      if (hit) return hit;
      const r = await fetch(e.request);
      if (r.ok) c.put(e.request, r.clone());
      return r;
    }));
    return;
  }
  // latest.json, la app y Open-Meteo -> primero red, si falla la copia guardada
  if (url.origin === location.origin || url.hostname.endsWith("open-meteo.com")) {
    e.respondWith(fetch(e.request).then(r => {
      if (r.ok && url.origin === location.origin) { const copia = r.clone(); caches.open(VERSION).then(c => c.put(e.request, copia)); }
      return r;
    }).catch(() => caches.match(e.request, {ignoreSearch: url.pathname.endsWith("latest.json")})));
  }
});
