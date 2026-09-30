/* Бурсач: офлайн-кэш. Версия меняется при каждой сборке. */
var V = "bursach-48510b6a", CORE = ["./", "index.html", "app.js?v=48510b6a", "manifest.json", "icon-180.png", "icon-192.png", "icon-512.png"];
self.addEventListener("install", function (e) { self.skipWaiting(); e.waitUntil(caches.open(V).then(function (c) { return c.addAll(CORE); })); });
self.addEventListener("activate", function (e) { e.waitUntil(caches.keys().then(function (ks) { return Promise.all(ks.filter(function (k) { return k !== V; }).map(function (k) { return caches.delete(k); })); }).then(function () { return self.clients.claim(); })); });
self.addEventListener("fetch", function (e) {
  var r = e.request; if (r.method !== "GET") return;
  var u = new URL(r.url);
  if (r.mode === "navigate") {                       /* страница: сеть с таймаутом 3 с, иначе кэш */
    e.respondWith(new Promise(function (res) {
      var done = false, t = setTimeout(function () { caches.match("index.html").then(function (m) { if (m && !done) { done = true; res(m); } }); }, 3000);
      fetch(r).then(function (n) { var cp = n.clone(); caches.open(V).then(function (c) { c.put("index.html", cp); }); if (!done) { done = true; clearTimeout(t); res(n); } })
        .catch(function () { caches.match("index.html").then(function (m) { if (!done) { done = true; res(m || Response.error()); } }); });
    })); return;
  }
  if (u.origin === location.origin || /fonts\.(googleapis|gstatic)\.com$/.test(u.hostname)) {   /* остальное: кэш, затем сеть */
    e.respondWith(caches.match(r).then(function (m) { return m || fetch(r).then(function (n) { if (n.ok || n.type === "opaque") { var cp = n.clone(); caches.open(V).then(function (c) { c.put(r, cp); }); } return n; }); }));
  }
});
