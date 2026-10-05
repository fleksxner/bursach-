/* ======================================================================
   ГИДРОМЕТЕОРОЛОГИЧЕСКОЕ ОБЕСПЕЧЕНИЕ СУДОВОЖДЕНИЯ (ГМОС)
   ПР №1 — атмосферное давление, истинный ветер, температура и влажность
           (половинка листа А4; варианты выдают на занятии — числа вводятся с листа).
   ПР №3 — КР №3: глубина под килем по Таблицам приливов, приливо-отливное,
           ветровое и суммарное течение.
   Константы приборов (сертификат анероида, график МС-13) восстановлены по
   проверенным преподавателем работам и редактируются в форме.
   Психрометрические таблицы (1981) воспроизведены их же формулами:
   осн. таблица — станционный психрометр при 1000 гПа:
     eт = E(t′) − 0,797·(t − t′)·(1 + 0,00115·t′),  E — по Тетенсу (над водой);
   табл. 4 (аспирационный психрометр): Δe = (t − t′)·(0,7947 − 0,000662·P).
   Исходники: папка «гмос/исходники приложения».
   ====================================================================== */
var GMOS = (function () {
  "use strict";
  var DISC = App.discipline("gmos", "Гидрометеорологическое обеспечение судовождения", "ГМОС");
  var R = Math.PI / 180;
  function num(x) { var s = String(x == null ? "" : x).trim().replace(/[−–]/g, "-").replace(",", "."); return s === "" ? NaN : parseFloat(s); }
  function f(x, n) { return App.f(x, n === undefined ? 1 : n).replace(/^-/, "−"); }
  function fs(x, n) { if (!isFinite(x)) return "—"; var s = f(x, n); if (s.charAt(0) === "−" || /^0(,0+)?$/.test(s)) return s; return "+" + s; }
  function fp(x, n) { var t = f(x, n); return t.charAt(0) === "−" ? "(" + t + ")" : t; }
  function deg(x) { x = ((x % 360) + 360) % 360; var r = Math.round(x); return r === 0 ? 360 : r; }
  function norm(x) { return ((x % 360) + 360) % 360; }
  /* строки таблицы: числа через пробел / «;» / «|» / «/»; запятая — десятичный знак */
  function lines(txt) { return String(txt || "").split(/\n+/).map(function (s) { return s.trim(); }).filter(function (s) { return s && !/^#/.test(s); }); }
  function toks(s) { return s.replace(/[−–]/g, "-").split(/[\s;|\/\t]+/).filter(Boolean); }
  function isNum(t) { return /^[+\-]?\d+(?:[.,]\d+)?$/.test(t); }
  /* румб → градусы */
  var RUMB = { N: 0, NNE: 22.5, NE: 45, ENE: 67.5, E: 90, ESE: 112.5, SE: 135, SSE: 157.5, S: 180, SSW: 202.5, SW: 225, WSW: 247.5, W: 270, WNW: 292.5, NW: 315, NNW: 337.5 };
  function dirOf(t) {
    if (isNum(t)) return num(t);
    var s = String(t).toUpperCase().replace(/Е/g, "E").replace(/Н/g, "N");
    if (/^[СЮВЗ]+$/.test(s)) s = s.replace(/С/g, "N").replace(/Ю/g, "S").replace(/В/g, "E").replace(/З/g, "W");
    return RUMB.hasOwnProperty(s) ? RUMB[s] : NaN;
  }
  /* таблица «x y» → линейная интерполяция (с экстраполяцией по крайнему отрезку) */
  function ptsOf(txt) { var p = []; lines(txt).forEach(function (l) { var t = toks(l).filter(isNum); if (t.length >= 2) p.push([num(t[0]), num(t[1])]); }); p.sort(function (a, b) { return a[0] - b[0]; }); return p; }
  function interp(p, x) {
    if (!p.length) return { v: NaN };
    if (p.length === 1) return { v: p[0][1], a: p[0], b: p[0], out: x !== p[0][0] };
    var i = 1; while (i < p.length - 1 && x > p[i][0]) i++;
    var a = p[i - 1], b = p[i]; return { v: a[1] + (b[1] - a[1]) * (x - a[0]) / (b[0] - a[0]), a: a, b: b, out: x < p[0][0] || x > p[p.length - 1][0] };
  }
  /* диапазоны «от до поправка» */
  function rngOf(txt) { var r = []; lines(txt).forEach(function (l) { var t = toks(l).filter(isNum); if (t.length >= 3) r.push([num(t[0]), num(t[1]), num(t[2])]); }); return r; }
  function rngCorr(r, x) { for (var i = 0; i < r.length; i++) if (x >= r[i][0] - 1e-9 && x <= r[i][1] + 1e-9) return { v: r[i][2], r: r[i] }; return { v: 0, r: null }; }

  /* ---------- психрометрия ---------- */
  function Ew(t) { return 6.1078 * Math.pow(10, 7.5 * t / (237.3 + t)); }          /* над водой, гПа */
  function Ei(t) { return 6.1078 * Math.pow(10, 9.5 * t / (265.5 + t)); }          /* надо льдом */
  var PSY = { Aw: 0.797, Ai: 0.702, Aaw: 0.7947, Bw: 0.000662, Aai: 0.702, Bi: 0.000585 };
  function eTab(t, tw, ice) { return (ice ? Ei(tw) : Ew(tw)) - (ice ? PSY.Ai : PSY.Aw) * (t - tw) * (1 + 0.00115 * tw); }
  function dE(dt, P, ice) { return dt * (ice ? PSY.Aai - PSY.Bi * P : PSY.Aaw - PSY.Bw * P); }
  function solveT(fn, y, lo, hi) { for (var k = 0; k < 80; k++) { var m = (lo + hi) / 2; if (fn(m) < y) lo = m; else hi = m; } return (lo + hi) / 2; }

  /* ---------- общее для рисунков ---------- */
  function arrow(c, x0, y0, x1, y1, col, w, head) {
    c.strokeStyle = col; c.fillStyle = col; c.lineWidth = w || 3; c.beginPath(); c.moveTo(x0, y0); c.lineTo(x1, y1); c.stroke();
    var a = Math.atan2(y1 - y0, x1 - x0), h = head || 18; c.beginPath(); c.moveTo(x1, y1);
    c.lineTo(x1 - h * Math.cos(a - 0.33), y1 - h * Math.sin(a - 0.33)); c.lineTo(x1 - h * Math.cos(a + 0.33), y1 - h * Math.sin(a + 0.33)); c.closePath(); c.fill();
  }
  function dimV(c, x, ya, yb, txt, col, FONT, left) {   /* вертикальный размер со стрелками */
    if (!isFinite(ya) || !isFinite(yb) || Math.abs(yb - ya) < 2) return;
    c.strokeStyle = col; c.fillStyle = col; c.lineWidth = 2; c.beginPath(); c.moveTo(x, ya); c.lineTo(x, yb); c.stroke();
    [[ya, yb], [yb, ya]].forEach(function (p) { var s = p[1] > p[0] ? 1 : -1; c.beginPath(); c.moveTo(x, p[0]); c.lineTo(x - 6, p[0] + 13 * s); c.lineTo(x + 6, p[0] + 13 * s); c.closePath(); c.fill(); });
    c.font = "600 21px " + FONT; c.textAlign = left ? "right" : "left"; c.fillText(txt, x + (left ? -10 : 10), (ya + yb) / 2 + 7);
  }
  function compass(c, cx, cy, rr, FONT, step, lab) {  /* круг с градусными делениями */
    c.strokeStyle = "#333"; c.lineWidth = 2; c.beginPath(); c.arc(cx, cy, rr, 0, 2 * Math.PI); c.stroke();
    c.font = "17px " + FONT; c.fillStyle = "#333"; c.textAlign = "center";
    for (var a = 0; a < 360; a += step) {
      var big = a % (lab || 30) === 0, s = Math.sin(a * R), co = Math.cos(a * R);
      c.lineWidth = big ? 2 : 1; c.beginPath(); c.moveTo(cx + rr * s, cy - rr * co); c.lineTo(cx + (rr - (big ? 16 : 8)) * s, cy - (rr - (big ? 16 : 8)) * co); c.stroke();
      if (big) c.fillText(String(a), cx + (rr + 22) * s, cy - (rr + 22) * co + 6);
    }
  }
  function P2(cx, cy, a, l) { return [cx + l * Math.sin(a * R), cy - l * Math.cos(a * R)]; }

  /* ======================================================================
     ПР №1
     ====================================================================== */
  var DEF1 = {
    variant: "1", vno: "2", H: "15",
    pr: "1 18,6 770,1\n4 8,4 749,7\n8 -1,2 742,9",
    wd: "1 235 16 270 8344 6992\n4 223 13 170 3480 2350\n5 265 16 310 4715 3470\n9 285 12 NW 12\n10 224 17 S 17",
    hm: "1 24,7 17,6 1000\n8 0,4 -3,2 1020\n5 4,8 2,5 1012",
    c: "-0,03", dpd: "-1,1",
    shk: "730 +2,1\n740 +1,4\n750 0,0\n760 -2,1\n770 -3,6\n780 -1,9\n790 -2,0",
    ma: "0,52", mb: "0,86", mt: "100", kn: "0,514", meth: "СМО СМО СМО граф аналит",
    tdry: "0,1 4,9 +0,6\n5,0 9,0 +0,3\n9,1 14,1 +0,4\n14,2 17,3 +0,2\n17,4 21,6 +0,1\n21,7 30,2 +0,2",
    twet: "0,0 4,8 +0,5\n4,9 9,3 +0,2\n9,4 13,9 -0,3\n14,0 16,9 -0,5\n17,0 22,1 -0,2\n22,2 29,2 +0,1",
    bat: "w"
  };
  function parseRows(txt, nNum) {   /* строки с № в начале (если чисел больше, чем нужно) */
    var out = [], k = 0;
    lines(txt).forEach(function (l) {
      var t = toks(l), ns = t.filter(isNum); k++;
      if (ns.length < nNum) return;
      var no = ns.length > nNum ? ns[0] : String(k), v = ns.slice(ns.length > nNum ? 1 : 0).map(num);
      out.push({ no: no, v: v, raw: l });
    });
    return out;
  }
  function pr1Pressure(d, st, warn) {
    var c = num(d.c), dpd = num(d.dpd), H = num(d.H), shk = ptsOf(d.shk), rows = parseRows(d.pr, 2), res = [];
    if (!shk.length) warn.push("Не задана таблица шкаловой поправки анероида — ΔPшк принята 0.");
    rows.forEach(function (r) {
      var t = r.v[0], P = r.v[1], dPt = Math.round(c * t * 1000) / 1000, ip = interp(shk, P), dPs = shk.length ? Math.round(ip.v * 100) / 100 : 0;
      if (ip.out) warn.push("Давление " + f(P, 1) + " мм вне таблицы шкаловой поправки — поправка экстраполирована.");
      var sum = Math.round((dPt + dpd + dPs) * 1000) / 1000, Pi = Math.round((P + sum) * 100) / 100;
      var h = Math.round(8000 / Pi * (1 + 0.004 * t) * 100) / 100, Hh = Math.round(H / h * 100) / 100, P0 = Math.round((Pi + Hh) * 100) / 100;
      res.push({ no: r.no, t: t, P: P, dPt: dPt, dpd: dpd, dPs: dPs, ip: ip, sum: sum, Pi: Pi, h: h, Hh: Hh, P0: P0 });
    });
    var L = ["ΔPt = c·t, c = " + f(c, 2) + " мм/°С (температурный коэффициент анероида);  ΔPдоб = " + f(dpd, 1) + " мм (добавочная поправка);  ΔPшк — по таблице шкаловых поправок (линейная интерполяция)",
      "ΣP = ΔPt + ΔPдоб + ΔPшк;  Pиспр = P + ΣP",
      "Барическая ступень h = (8000 / Pиспр)·(1 + 0,004·t), м/мм рт. ст.;  Pморя = Pиспр + H/h,  H = " + f(H, 0) + " м"];
    res.forEach(function (x) {
      L.push("№ " + x.no + ":  ΔPt = " + f(c, 2) + "·" + fp(x.t, 1) + " = " + fs(x.dPt, 3) + ";  ΔPшк(" + f(x.P, 1) + ") = " + fs(x.dPs, 2) +
        (x.ip.a && x.ip.b && x.ip.a !== x.ip.b ? " [между " + f(x.ip.a[0], 0) + " → " + fs(x.ip.a[1], 1) + " и " + f(x.ip.b[0], 0) + " → " + fs(x.ip.b[1], 1) + "]" : "") +
        ";  ΣP = " + fs(x.dPt, 3) + " " + (x.dpd < 0 ? "− " : "+ ") + f(Math.abs(x.dpd), 1) + " " + (x.dPs < 0 ? "− " : "+ ") + f(Math.abs(x.dPs), 2) + " = " + fs(x.sum, 3) +
        ";  Pиспр = " + f(x.P, 1) + " " + (x.sum < 0 ? "− " : "+ ") + f(Math.abs(x.sum), 3) + " = " + f(x.Pi, 2) + " мм");
      L.push("        h = 8000/" + f(x.Pi, 2) + "·(1 + 0,004·" + fp(x.t, 1) + ") = " + f(x.h, 2) + " м/мм;  H/h = " + f(H, 0) + "/" + f(x.h, 2) + " = " + f(x.Hh, 2) +
        ";  Pморя = " + f(x.Pi, 2) + " + " + f(x.Hh, 2) + " = " + f(x.P0, 2) + " мм рт. ст. (≈ " + f(x.P0 * 1.33322, 1) + " гПа)");
    });
    st.push({ no: 1, title: "Определение атмосферного давления и приведение его к уровню моря (" + f(H, 0) + " м)", lines: [], after: L,
      tables: [{ headers: ["№", "tсух.испр, °С", "P, мм рт.ст.", "ΔPt", "ΔPдоб", "ΔPшк", "ΣP", "Pиспр", "h", "H/h", "Pморя"],
        rows: res.map(function (x) { return [x.no, f(x.t, 1), f(x.P, 1), fs(x.dPt, 3), fs(x.dpd, 1), fs(x.dPs, 2), fs(x.sum, 3), f(x.Pi, 2), f(x.h, 2), f(x.Hh, 2), f(x.P0, 2)]; }),
        widths: [0.9, 1.6, 1.7, 1.5, 1.3, 1.4, 1.5, 1.8, 1.3, 1.2, 1.8], size: 9 }],
      answer: res.length ? "Pморя: " + res.map(function (x) { return "№ " + x.no + " — " + f(x.P0, 2); }).join("; ") + " мм рт. ст." : "" });
    if (!res.length) warn.push("Задача 1: введите строки «№ t P» (например «1 18,6 770,1»).");
    return res;
  }
  function trueWind(K, Vs, A, W) {
    var x = W * Math.sin(A * R) - Vs * Math.sin(K * R), y = W * Math.cos(A * R) - Vs * Math.cos(K * R);
    var U = Math.sqrt(x * x + y * y), Ku = norm(Math.atan2(x, y) / R);
    var q = norm(A - K); if (q > 180) q -= 360;                         /* КУ кажущегося ветра: + правый борт */
    var Ua = Math.sqrt(Vs * Vs - 2 * W * Vs * Math.cos(q * R) + W * W);
    var cqi = Ua > 1e-9 ? (W * Math.cos(q * R) - Vs) / Ua : 1, qi = Math.acos(Math.max(-1, Math.min(1, cqi))) / R;
    var sg = Math.sin(q * R) < 0 ? -1 : 1;
    return { U: U, K: Ku, q: q, Ua: Ua, cqi: cqi, qi: qi, sg: sg };
  }
  function pr1Wind(d, st, warn) {
    var a = num(d.ma), b = num(d.mb), tt = num(d.mt) || 100, kn = num(d.kn) || 0.514, M = String(d.meth || "").split(/[\s,;]+/).filter(Boolean), res = [], k = 0;
    lines(d.wd).forEach(function (l) {
      var t = toks(l), meth = null; k++;
      t = t.filter(function (x) { if (/^(смо|граф|анал)/i.test(x)) { meth = /^смо/i.test(x) ? "СМО" : /^граф/i.test(x) ? "граф" : "аналит"; return false; } return true; });
      var dirI = -1; t.forEach(function (x, i) { if (!isNum(x) && dirI < 0) dirI = i; });
      var r = { raw: l }, v;
      if (dirI >= 0) {           /* направление румбом: [№] K V румб (II I | АРИ) */
        var pre = t.slice(0, dirI).filter(isNum).map(num), post = t.slice(dirI + 1).filter(isNum).map(num);
        r.no = pre.length >= 3 ? String(pre[0]) : String(k); pre = pre.slice(-2);
        r.K = pre[0]; r.Vk = pre[1]; r.A = dirOf(t[dirI]); r.Atxt = t[dirI].toUpperCase().replace(/Е/g, "E");
        if (post.length >= 2) { r.II = post[0]; r.I = post[1]; } else r.W0 = post[0];
      } else {                   /* всё числами */
        v = t.filter(isNum).map(num);
        if (v.length >= 6) { r.no = String(v[0]); r.K = v[1]; r.Vk = v[2]; r.A = v[3]; r.II = v[4]; r.I = v[5]; }
        else if (v.length === 5) {
          if (v[3] > 60 && v[4] > 60) { r.no = String(k); r.K = v[0]; r.Vk = v[1]; r.A = v[2]; r.II = v[3]; r.I = v[4]; }
          else { r.no = String(v[0]); r.K = v[1]; r.Vk = v[2]; r.A = v[3]; r.W0 = v[4]; }
        } else if (v.length === 4) { r.no = String(k); r.K = v[0]; r.Vk = v[1]; r.A = v[2]; r.W0 = v[3]; }
        r.Atxt = isFinite(r.A) ? f(r.A, 0) : "?";
      }
      if (![r.K, r.Vk, r.A].every(isFinite)) { warn.push("Задача 2: не разобрана строка «" + l + "»."); return; }
      if (isFinite(r.II) && isFinite(r.I)) {
        var dn = r.II - r.I; if (dn < 0) dn += 10000; r.n = Math.round(dn / tt * 100) / 100; r.W = Math.round((a + b * r.n) * 10) / 10; r.mc = true;
      } else if (isFinite(r.W0)) { r.W = r.W0; r.mc = false; }
      else { warn.push("Задача 2: в строке «" + l + "» нет отсчётов анемометра."); return; }
      r.meth = meth || M[res.length] || (res.length < 3 ? "СМО" : res.length === 3 ? "граф" : "аналит");
      r.Vs = Math.round(r.Vk * kn * 10) / 10;
      r.tw = trueWind(r.K, r.Vs, r.A, r.W);
      res.push(r);
    });
    var L = ["Скорость судна в м/с: V = Vуз·" + f(kn, 3) + " (1 уз = 1852 м/ч).  МС-13: n = (II − I)/" + f(tt, 0) + " с (дел/с), Wк = " + f(a, 2) + " + " + f(b, 2) + "·n (график сертификата анемометра)",
      "Аналитически: q = Aк − ИК (КУ кажущегося ветра), Vи = √(V² − 2·Wк·V·cos q + Wк²), cos qи = (Wк·cos q − V)/Vи, Kи = ИК ± qи (знак как у q: + правый борт)"];
    res.forEach(function (r) {
      var w = r.tw, s = "№ " + r.no + " (" + r.meth + "):  V = " + f(r.Vk, 0) + "·" + f(kn, 3) + " = " + f(r.Vs, 1) + " м/с;  " +
        (r.mc ? "n = (" + r.II + " − " + r.I + ")/" + f(tt, 0) + " = " + f(r.n, 2) + " дел/с → Wк = " + f(a, 2) + " + " + f(b, 2) + "·" + f(r.n, 2) + " = " + f(r.W, 1) + " м/с" : "Wк = " + f(r.W, 1) + " м/с (АРИ-49)") +
        ";  кажущийся ветер " + r.Atxt + (isNaN(+r.Atxt) ? " (" + f(r.A, 1).replace(/,0$/, "") + "°)" : "°");
      L.push(s);
      L.push("        q = " + f(r.A, 1).replace(/,0$/, "") + " − " + f(r.K, 0) + " = " + fs(w.q, 1).replace(/,0$/, "") + "°;  Vи = √(" + f(r.Vs, 1) + "² − 2·" + f(r.W, 1) + "·" + f(r.Vs, 1) + "·cos " + f(Math.abs(w.q), 1).replace(/,0$/, "") + "° + " + f(r.W, 1) + "²) = " + f(w.Ua, 1) +
        " м/с;  cos qи = (" + f(r.W, 1) + "·cos " + f(Math.abs(w.q), 1).replace(/,0$/, "") + "° − " + f(r.Vs, 1) + ")/" + f(w.Ua, 2) + " = " + f(w.cqi, 3) + " → qи = " + f(w.qi, 1) + "°;  Kи = " + f(r.K, 0) + (w.sg < 0 ? " − " : " + ") + f(w.qi, 1) + " = " + deg(r.K + w.sg * w.qi) + "°");
      if (r.meth === "СМО") L.push("        Круг СМО: ИК " + f(r.K, 0) + "° к индексу, от центра к индексу ОК = " + f(r.Vs, 1) + " м/с; к индексу " + f(r.A, 1).replace(/,0$/, "") + "°, ОВ = " + f(r.W, 1) +
        " м/с; вращать круг, пока В и К не лягут на линию, параллельную диаметру через индекс (В ниже К): ВК = " + f(w.U, 1) + " м/с, у индекса " + deg(w.K) + "°.");
      if (r.meth === "граф") L.push("        Графически (рис.): ОК = " + f(r.Vs, 1) + " м/с по ИК " + f(r.K, 0) + "°, ОВ = " + f(r.W, 1) + " м/с по " + f(r.A, 1).replace(/,0$/, "") + "°; ОД ∥ КВ, ОД = КВ = " + f(w.U, 1) + " м/с, направление ОД = " + deg(w.K) + "°.");
    });
    st.push({ no: 2, title: "Определение истинного ветра на судне (Kи, Vи)", lines: [], after: L, figs: res.length ? ["wind"] : [],
      tables: [{ headers: ["№", "Курс, °", "V, уз", "Напр. каж. ветра", "II отсчёт", "I отсчёт", "Разность, дел/с", "АРИ-49 / Wк, м/с", "Kи, °", "Vи, м/с", "Способ"],
        rows: res.map(function (r) { return [r.no, f(r.K, 0), f(r.Vk, 0), r.Atxt, r.mc ? String(r.II) : "—", r.mc ? String(r.I) : "—", r.mc ? f(r.n, 2) : "—", f(r.W, 1), String(deg(r.tw.K)), f(r.tw.U, 1), r.meth]; }),
        widths: [0.8, 1.2, 1.1, 1.7, 1.5, 1.5, 1.6, 1.7, 1.2, 1.3, 1.4], size: 9 }],
      answer: res.length ? res.map(function (r) { return "№ " + r.no + ": " + deg(r.tw.K) + "° — " + f(r.tw.U, 1) + " м/с"; }).join("; ") : "" });
    if (!res.length) warn.push("Задача 2: введите строки «№ курс скорость направление II I» или «№ курс скорость направление АРИ».");
    return res;
  }
  function pr1Hum(d, st, warn) {
    var rd = rngOf(d.tdry), rw = rngOf(d.twet), ice0 = d.bat === "i", rows = parseRows(d.hm, 3), res = [], nz = [];
    rows.forEach(function (r) {
      var t = r.v[0], tw = r.v[1], P = r.v[2], cd = rngCorr(rd, t), cw = rngCorr(rw, tw);
      if (!cd.r || !cw.r) nz.push("№ " + r.no + ": " + (!cd.r ? "tсух = " + f(t, 1) : "") + (!cd.r && !cw.r ? ", " : "") + (!cw.r ? "tсмоч = " + f(tw, 1) : "") + " — в таблице поправок нет такого диапазона, поправка 0 (как принято в работах)");
      var ti = Math.round((t + cd.v) * 10) / 10, twi = Math.round((tw + cw.v) * 10) / 10, dt = Math.round((ti - twi) * 10) / 10, ice = ice0 && twi < 0;
      var et = Math.round(eTab(ti, twi, ice) * 10) / 10, de = Math.round(dE(dt, P, ice) * 100) / 100, e = Math.round((et + de) * 100) / 100;
      var tp = solveT(function (x) { return eTab(ti, x, ice); }, e, -60, ti), E = Ew(ti);
      var td = solveT(Ew, e, -80, 60), fr = 100 * e / E, dd = E - e, aa = 220 * e / (273 + ti);
      res.push({ no: r.no, t: t, tw: tw, P: P, cd: cd.v, cw: cw.v, ti: ti, twi: twi, dt: dt, et: et, de: de, e: e, tp: tp, td: td, f: fr, d: dd, E: E, a: aa, ice: ice });
    });
    var L = ["Поправки термометров — из свидетельства о поверке (диапазоны в форме);  Δt = tсух.испр − tсмоч.испр",
      "eт — осн. таблица Психрометрических таблиц по tсух.испр и tсмоч.испр;  Δe — табл. 4 (аспирационный психрометр) по P и Δt;  e = eт + Δe",
      "По e в столбце tсух.испр: t′ (строка таблицы с этим e), точка росы td, относительная влажность f = e/E·100 %, дефицит d = E − e;  абсолютная влажность a = 220·e/(273 + t)"].concat(nz);
    res.forEach(function (x) {
      L.push("№ " + x.no + ":  tсух = " + f(x.t, 1) + " " + (x.cd < 0 ? "− " : "+ ") + f(Math.abs(x.cd), 1) + " = " + f(x.ti, 1) + ";  tсмоч = " + f(x.tw, 1) + " " + (x.cw < 0 ? "− " : "+ ") + f(Math.abs(x.cw), 1) + " = " + f(x.twi, 1) +
        ";  Δt = " + f(x.dt, 1) + ";  eт = " + f(x.et, 1) + " гПа;  Δe(P = " + f(x.P, 0) + ", Δt = " + f(x.dt, 1) + ") = " + fs(x.de, 2) + ";  e = " + f(x.et, 1) + " + " + f(x.de, 2) + " = " + f(x.e, 2) + " гПа" + (x.ice ? " (батист — лёд)" : ""));
      L.push("        E(" + f(x.ti, 1) + ") = " + f(x.E, 1) + " гПа;  t′ = " + f(x.tp, 1) + ";  td = " + f(x.td, 1) + " °С;  f = " + f(x.e, 2) + "/" + f(x.E, 1) + "·100 = " + f(x.f, 0) + " %;  d = " + f(x.E, 1) + " − " + f(x.e, 2) + " = " + f(x.d, 1) + " гПа;  a = 220·" + f(x.e, 2) + "/" + f(273 + x.ti, 1) + " = " + f(x.a, 1) + " г/м³");
    });
    st.push({ no: 3, title: "Измерение температуры и влажности на судне (по Психрометрическим таблицам)", lines: [], after: L,
      tables: [{ headers: ["№", "tсух, °С", "tсмоч, °С", "P, гПа", "Δtсух", "Δtсмоч", "tсух испр", "tсмоч испр", "Δt", "eт", "Δe", "e"],
        rows: res.map(function (x) { return [x.no, f(x.t, 1), f(x.tw, 1), f(x.P, 0), fs(x.cd, 1), fs(x.cw, 1), f(x.ti, 1), f(x.twi, 1), f(x.dt, 1), f(x.et, 1), f(x.de, 2), f(x.e, 2)]; }),
        widths: [0.8, 1.3, 1.3, 1.2, 1.2, 1.3, 1.4, 1.5, 1.0, 1.2, 1.2, 1.3], size: 9 },
      { headers: ["№", "t′", "td", "e", "f, %", "d"], rows: res.map(function (x) { return [x.no, f(x.tp, 1), f(x.td, 1), f(x.e, 2), f(x.f, 0), f(x.d, 1)]; }), widths: [1, 2.5, 2.5, 2.5, 2.5, 2.5], size: 9.5 }],
      answer: res.length ? res.map(function (x) { return "№ " + x.no + ": e = " + f(x.e, 2) + " гПа, td = " + f(x.td, 1) + "°, f = " + f(x.f, 0) + " %, d = " + f(x.d, 1) + " гПа"; }).join("; ") : "" });
    if (!res.length) warn.push("Задача 3: введите строки «№ tсух tсмоч P».");
    return res;
  }
  function drawWind(cv, LK, d, r) {
    var W = 1000, H = 1130; cv.width = W; cv.height = H; var c = cv.getContext("2d"), FONT = LK.font;
    c.fillStyle = "#fff"; c.fillRect(0, 0, W, H);
    var rows = r.wind || [], x = rows.filter(function (z) { return z.meth === "граф"; })[0] || rows[0];
    if (!x) return;
    var cx = W / 2, cy = 460, RR = 380, mx = Math.max(x.Vs, x.W, x.tw.U, 1), sc = (RR - 40) / mx;
    compass(c, cx, cy, RR, FONT, 5, 30);
    c.strokeStyle = "#ddd"; c.lineWidth = 1; for (var k = 1; k <= Math.floor(mx); k++) { c.beginPath(); c.arc(cx, cy, k * sc, 0, 2 * Math.PI); c.stroke(); }
    var K = P2(cx, cy, x.K, x.Vs * sc), B = P2(cx, cy, x.A, x.W * sc), D = [cx + B[0] - K[0], cy + B[1] - K[1]];
    c.setLineDash([8, 7]); c.strokeStyle = "#888"; c.lineWidth = 2; c.beginPath(); c.moveTo(K[0], K[1]); c.lineTo(B[0], B[1]); c.lineTo(D[0], D[1]); c.stroke(); c.setLineDash([]);
    arrow(c, cx, cy, K[0], K[1], "#1d4e89", 4); arrow(c, cx, cy, B[0], B[1], "#2a9d8f", 4); arrow(c, cx, cy, D[0], D[1], "#c1121f", 5);
    c.font = "600 24px " + FONT; c.textAlign = "center";
    c.fillStyle = "#1d4e89"; c.fillText("К", K[0] + 22 * Math.sin(x.K * R), K[1] - 22 * Math.cos(x.K * R) + 8);
    c.fillStyle = "#2a9d8f"; c.fillText("В", B[0] + 22 * Math.sin(x.A * R), B[1] - 22 * Math.cos(x.A * R) + 8);
    c.fillStyle = "#c1121f"; c.fillText("Д", D[0] + 24 * Math.sin(x.tw.K * R), D[1] - 24 * Math.cos(x.tw.K * R) + 8);
    c.fillStyle = "#111"; c.fillText("О", cx - 18, cy + 26);
    c.font = "20px " + FONT; c.textAlign = "left"; c.fillStyle = "#111";
    var tx = 24, ty = 30;
    ["Строка № " + x.no + " (" + x.meth + "), масштаб: 1 круг = 1 м/с", "ОК — курсовой ветер: ИК " + f(x.K, 0) + "°, V = " + f(x.Vs, 1) + " м/с", "ОВ — кажущийся ветер: " + f(x.A, 1).replace(/,0$/, "") + "°, Wк = " + f(x.W, 1) + " м/с",
      "ОД ∥ КВ — истинный ветер: Kи = " + deg(x.tw.K) + "°, Vи = " + f(x.tw.U, 1) + " м/с"].forEach(function (s, i) {
      c.fillStyle = ["#111", "#1d4e89", "#2a9d8f", "#c1121f"][i]; c.fillText(s, tx, H - 128 + i * 30);
    });
  }
  App.taskWork({
    disc: "gmos", discName: "Гидрометеорологическое обеспечение судовождения", id: "gmos-pr1", order: 1, key: "gmos-pr1", file: "ГМОС_ПР1",
    no: "ПР 1", short: "ПР №1", eyebrow: "ГМОС · Практическая работа №1", titleKind: "Практическая работа",
    title: "Атмосферное давление, истинный ветер, температура и влажность на судне",
    titleTopic: "Определение атмосферного давления, истинного ветра, температуры и влажности воздуха на судне",
    desc: "Половинка листа: давление → уровень моря (15 м), истинный ветер (СМО, график, аналитически), влажность по Психрометрическим таблицам. Числа вводятся с листа варианта.",
    variants: 1, variantData: function () { return Object.assign({}, DEF1); },
    varText: function (d) { return "Вариант " + (d.vno || "—"); },
    hint: "Перепиши числа со своего листа: каждая строка таблицы — одна строка поля, числа через пробел, № строки — первым. По умолчанию стоит вариант № 2 (проверенный преподавателем). Данные приборов внизу — как давали на занятиях; если на доске были другие, поправь.",
    sections: [
      { title: "Лист варианта", fields: [["vno", "Вариант №"], ["H", "Высота рубки H", { num: true, unit: "м" }]] },
      { title: "1. Давление", note: "Строка: «№ tсух.испр P(мм рт. ст.)», например «1 18,6 770,1».", fields: [["pr", "№  t, °С  P, мм", { area: true, rows: 3 }]] },
      { title: "2. Истинный ветер", note: "МС-13: «№ курс скорость направление II I» (например «4 168 13 280 6164 5453»). АРИ-49: «№ курс скорость направление отсчёт» (например «6 138 15 SE 7»). Направление — градусы или румб (N, NE, SE…). Способ — по порядку строк (поле внизу) или словом в конце строки: СМО / граф / аналит.",
        fields: [["wd", "Строки задачи 2", { area: true, rows: 5 }]] },
      { title: "3. Температура и влажность", note: "Строка: «№ tсух tсмоч P(гПа)», например «2 21,1 17,2 1000».", fields: [["hm", "№  tсух  tсмоч  P", { area: true, rows: 3 }], ["bat", "Батист при t′ < 0", { select: [["w", "вода (как в работах)"], ["i", "лёд"]] }]] },
      { title: "Данные приборов (с занятий)", note: "Восстановлены по проверенным работам; поправь, если на доске были другие.", fields: [
        ["c", "Темп. коэфф. анероида c", { num: true, unit: "мм/°С" }], ["dpd", "Добавочная поправка ΔPдоб", { num: true, unit: "мм" }],
        ["shk", "Шкаловая поправка: «P ΔPшк»", { area: true, rows: 7 }],
        ["ma", "МС-13: Wк = a + b·n, a", { num: true }], ["mb", "b", { num: true }], ["mt", "Время замера", { num: true, unit: "с" }], ["kn", "1 уз =", { num: true, unit: "м/с" }],
        ["meth", "Способы по строкам задачи 2", { wide: true }],
        ["tdry", "Поправки сухого термометра «от до поправка»", { area: true, rows: 6 }], ["twet", "Поправки смоченного термометра «от до поправка»", { area: true, rows: 6 }]] }],
    goal: "Цель: определить исправленное атмосферное давление и привести его к уровню моря, рассчитать направление и скорость истинного ветра, определить характеристики влажности воздуха по Психрометрическим таблицам.",
    figs: { wind: { caption: "Графическое определение истинного ветра", draw: drawWind } },
    solve: function (d) {
      var st = [], warn = [], p = pr1Pressure(d, st, warn), w = pr1Wind(d, st, warn), h = pr1Hum(d, st, warn);
      return { steps: st, warn: warn, press: p, wind: w, hum: h };
    }
  });

  return { num: num, f: f, fs: fs, deg: deg, norm: norm, lines: lines, toks: toks, isNum: isNum, arrow: arrow, dimV: dimV, compass: compass, P2: P2, R: R,
    eTab: eTab, dE: dE, Ew: Ew, trueWind: trueWind };
})();

/* ===== ГМОС · ПР №3 (КР №3): глубина под килем, приливо-отливное, ветровое и суммарное течение ===== */
(function () {
  "use strict";
  var G = GMOS, num = G.num, f = G.f, fs = G.fs, deg = G.deg, norm = G.norm, lines = G.lines, toks = G.toks, isNum = G.isNum, R = G.R;
  function pad(n) { return (n < 10 ? "0" : "") + n; }
  function hm(s) { var m = String(s == null ? "" : s).match(/(\d{1,2})\s*[:.,чh ]\s*(\d{1,2})/); if (m) return +m[1] + +m[2] / 60; var n = String(s || "").match(/^\s*(\d{1,2})\s*$/); return n ? +n[1] : NaN; }
  function T(h) { if (!isFinite(h)) return "—"; var d = Math.floor(h / 24 + 1e-9), x = h - d * 24, H = Math.floor(x + 1e-9), M = Math.round((x - H) * 60); if (M === 60) { H++; M = 0; } if (H === 24) { H = 0; d++; } return pad(H) + ":" + pad(M) + (d < 0 ? " (пред.)" : d > 0 ? " (след.)" : ""); }
  function dur(h) { var H = Math.floor(Math.abs(h) + 1e-9), M = Math.round((Math.abs(h) - H) * 60); if (M === 60) { H++; M = 0; } return (h < 0 ? "−" : "") + H + " ч " + pad(M) + " мин"; }
  function dT(s) {                   /* поправка времени «±ч.мм» (или «±мм») → часы */
    var t = String(s == null ? "" : s).trim().replace(/[−–]/g, "-"); if (!t || t === "-" || t === "—") return 0;
    var sg = /^-/.test(t) ? -1 : 1, m = t.match(/(\d{1,2})\s*[:.,]\s*(\d{1,2})/);
    if (m) return sg * (+m[1] + +m[2] / 60);
    var n = t.match(/(\d+)/); return n ? sg * +n[1] / 60 : 0;
  }
  function dTtxt(s) { var v = dT(s); if (!v) return "—"; var a = Math.abs(v), H = Math.floor(a + 1e-9), M = Math.round((a - H) * 60); if (M === 60) { H++; M = 0; } return (v < 0 ? "−" : "+") + H + "." + pad(M); }
  function hv(x) { var v = num(x); return isFinite(v) ? v : 0; }
  function parseEv(txt) {
    var ev = [];
    lines(txt).forEach(function (s) {
      var day = 0; if (/^(-|−|пред)/i.test(s)) day = -1; else if (/^(\+|след)/i.test(s)) day = 1;
      var m = s.replace(/^(-|−|\+|пред\S*|след\S*)\s*/i, "").match(/(\d{1,2})\s*[:.,]\s*(\d{2})\s+(-?\d+(?:[.,]\d+)?)/); if (!m) return;
      ev.push({ t: day * 24 + +m[1] + +m[2] / 60, h: num(m[3]) });
    });
    ev.sort(function (a, b) { return a.t - b.t; });
    ev.forEach(function (e, i) { var p = ev[i - 1], n = ev[i + 1]; e.pv = (p ? e.h > p.h : n ? e.h > n.h : true); });
    return ev;
  }
  function hAt(ev, t) {
    for (var i = 1; i < ev.length; i++) if (t >= ev[i - 1].t - 1e-9 && t <= ev[i].t + 1e-9) {
      var a = ev[i - 1], b = ev[i], k = (t - a.t) / (b.t - a.t); return a.h + (b.h - a.h) * (1 - Math.cos(Math.PI * k)) / 2;
    }
    return NaN;
  }

  var DEF3 = {
    variant: "1", vno: "21", hsel: "sz", P0: "1013", kP: "0,01",
    nm1: "б. Нагаева", T1: "12.00", dt1: "12.08", Hk1: "6,3", P1: "1034", Ts1: "6,6", no1: "86", co1: "", op1: "86, стр. 59",
    tp1: "", tm1: "", sp1: "", kp1: "", km1: "", sm1: "", ev1: "02.29 4,0\n08.48 1,9\n15.06 3,5\n20.31 2,1",
    nm2: "зал. Одян", T2: "14.00", dt2: "28.07", Hk2: "3,8", P2: "1028", Ts2: "6,2", no2: "91", co2: "", op2: "86, стр. 59",
    tp2: "-0.08", tm2: "-0.08", sp2: "+0,4", kp2: "+0,3", km2: "+0,2", sm2: "+0,1", ev2: "02.42 4,2\n09.01 1,5\n15.28 3,7\n21.04 1,9",
    cP: "пр. Акутан", cD: "07.10.2008", cT1: "8", cT2: "12", cN: "221,7", cC: "0,87", cPg: "89", cUseC: "1", cFl: "295",
    cH: "8 10 12", cRows: "222 59 295 31 295 14 115",
    phi: "43", hem: "N", kq: "330", vq: "1,8", kw: "170", w: "16", kc: "0,026"
  };
  function point(d, i, warn) {
    var p = { i: i, nm: d["nm" + i], T: hm(d["T" + i]), Ttxt: d["T" + i], dt: d["dt" + i], Hk: num(d["Hk" + i]), P: num(d["P" + i]), Ts: num(d["Ts" + i]),
      no: d["no" + i], co: d["co" + i], op: d["op" + i], tp: d["tp" + i], tm: d["tm" + i], sp: d["sp" + i], kp: d["kp" + i], km: d["km" + i], sm: d["sm" + i] };
    var sz = d.hsel !== "kv", P0 = num(d.P0) || 1013, kP = isFinite(num(d.kP)) ? num(d.kP) : 0.01;
    p.ev0 = parseEv(d["ev" + i]);
    p.cPV = hv(sz ? p.sp : p.kp); p.cMV = hv(sz ? p.sm : p.km); p.dPV = dT(p.tp); p.dMV = dT(p.tm);
    p.ev = p.ev0.map(function (e) { return { t: e.t + (e.pv ? p.dPV : p.dMV), h: Math.round((e.h + (e.pv ? p.cPV : p.cMV)) * 100) / 100, pv: e.pv, src: e }; });
    if (p.ev.length < 2) { warn.push("Пункт № " + i + ": введите не меньше двух ПВ/МВ основного пункта («чч.мм высота» по строкам)."); return p; }
    if (!isFinite(p.T)) { warn.push("Пункт № " + i + ": не задано время."); return p; }
    var seg = null; for (var k = 1; k < p.ev.length; k++) if (p.T >= p.ev[k - 1].t - 1e-9 && p.T <= p.ev[k].t + 1e-9) { seg = [p.ev[k - 1], p.ev[k]]; break; }
    if (!seg) { warn.push("Пункт № " + i + ": время " + T(p.T) + " вне введённых ПВ/МВ — добавьте воду предыдущих («-» в начале строки) или следующих («+») суток."); return p; }
    var a = seg[0], b = seg[1]; p.a = a; p.b = b;
    p.rise = b.h > a.h; p.pvW = a.pv ? a : b; p.mvW = a.pv ? b : a;
    p.Tpr = b.t - a.t; p.B = Math.round((p.pvW.h - p.mvW.h) * 100) / 100;
    p.near = Math.abs(p.T - a.t) <= Math.abs(b.t - p.T) ? a : b; p.dTb = Math.abs(p.T - p.near.t);
    p.dh = Math.round(p.B / 2 * (1 - Math.cos(Math.PI * p.dTb / p.Tpr)) * 100) / 100;
    p.hz = Math.round((p.near.pv ? p.near.h - p.dh : p.near.h + p.dh) * 100) / 100;
    p.dhp = Math.round((P0 - p.P) * kP * 100) / 100; p.P0 = P0; p.kP = kP;
    p.Hz = Math.round((p.Hk + p.hz + p.dhp - p.Ts) * 100) / 100;
    var mvs = p.ev.filter(function (e) { return !e.pv && e.t >= 0 && e.t < 24; }); if (!mvs.length) mvs = [p.mvW];
    p.mvMin = mvs.reduce(function (m, e) { return e.h < m.h ? e : m; }, mvs[0]);
    p.Hmin = Math.round((p.Hk + p.mvMin.h + p.dhp - p.Ts) * 100) / 100;
    if (p.Hmin < 0) warn.push("Пункт № " + i + ": Hмин = " + f(p.Hmin, 2) + " м < 0 — в малую воду " + T(p.mvMin.t) + " корабль коснётся грунта.");
    else if (p.Hz < 0) warn.push("Пункт № " + i + ": Hз < 0 — на заданный момент глубины под килем нет.");
    return p;
  }
  function evRows(ev) {          /* ПВ и МВ парами: [T ПВ, h ПВ, T МВ, h МВ] */
    var pv = ev.filter(function (e) { return e.pv; }), mv = ev.filter(function (e) { return !e.pv; }), n = Math.max(pv.length, mv.length), out = [];
    for (var k = 0; k < n; k++) out.push([pv[k] ? T(pv[k].t) : "—", pv[k] ? f(pv[k].h, 1) : "—", mv[k] ? T(mv[k].t) : "—", mv[k] ? f(mv[k].h, 1) : "—"]);
    return out;
  }
  /* ---------- приливо-отливное течение ---------- */
  function vec(v, k) { return [v * Math.sin(k * R), v * Math.cos(k * R)]; }
  function unvec(c) { return { v: Math.sqrt(c[0] * c[0] + c[1] * c[1]), k: norm(Math.atan2(c[0], c[1]) / R) }; }
  function current(d, warn) {
    var H = toks(String(d.cH || "")).filter(isNum).map(num), rows = [], N = num(d.cN), C = num(d.cC), useC = d.cUseC !== "0";
    lines(d.cRows).forEach(function (l) { var t = toks(l).filter(isNum).map(num); if (t.length >= 3) { var r = { N: t[0], val: [] }; for (var k = 1; k + 1 < t.length; k += 2) r.val.push({ v: t[k], k: t[k + 1] }); rows.push(r); } });
    var o = { H: H, rows: rows, N: N, C: C, useC: useC };
    if (!H.length || !rows.length) { warn.push("Течение: введите часы и строки таблицы «N V K V K …»."); return o; }
    rows.forEach(function (r) { if (r.val.length !== H.length) warn.push("Течение: в строке N = " + f(r.N, 1) + " пар «V K» — " + r.val.length + ", а часов — " + H.length + "."); });
    rows.sort(function (a, b) { return a.N - b.N; });
    var r1 = rows[0], r2 = rows[1] || null;
    if (rows.length > 2 && isFinite(N)) { for (var k = 1; k < rows.length; k++) if (N <= rows[k].N || k === rows.length - 1) { r1 = rows[k - 1]; r2 = rows[k]; break; } }
    o.r1 = r1; o.r2 = r2;
    var w = r2 && isFinite(N) && r2.N !== r1.N ? (N - r1.N) / (r2.N - r1.N) : 0;
    o.w = w;
    o.base = H.map(function (h, j) {
      var a = r1.val[j] || { v: 0, k: 0 }, b = r2 ? (r2.val[j] || a) : a, ca = vec(a.v, a.k), cb = vec(b.v, b.k);
      var u = unvec([ca[0] + (cb[0] - ca[0]) * w, ca[1] + (cb[1] - ca[1]) * w]);
      var k0 = u.v < 1e-9 ? a.k : u.k;
      return { h: h, v: Math.round(u.v * 10) / 10, k: deg(k0) };
    });
    o.base.forEach(function (x) { x.vc = useC && isFinite(C) ? Math.round(x.v * C * 10) / 10 : x.v; });
    /* по часам на заданный промежуток (нечётные часы — линейная интерполяция по времени) */
    var t1 = num(d.cT1), t2 = num(d.cT2); if (!isFinite(t1)) t1 = H[0]; if (!isFinite(t2)) t2 = H[H.length - 1];
    o.t1 = t1; o.t2 = t2; o.hours = [];
    for (var hh = Math.ceil(t1); hh <= Math.floor(t2) + 1e-9; hh++) {
      var ex = o.base.filter(function (x) { return x.h === hh; })[0];
      if (ex) { o.hours.push({ h: hh, v: ex.vc, k: ex.k, tab: true }); continue; }
      var lo = null, hi = null; o.base.forEach(function (x) { if (x.h < hh && (!lo || x.h > lo.h)) lo = x; if (x.h > hh && (!hi || x.h < hi.h)) hi = x; });
      if (!lo || !hi) { warn.push("Течение: час " + hh + " вне выписанных часов таблицы."); continue; }
      var q = (hh - lo.h) / (hi.h - lo.h), ca2 = vec(lo.vc, lo.k), cb2 = vec(hi.vc, hi.k), u2 = unvec([ca2[0] + (cb2[0] - ca2[0]) * q, ca2[1] + (cb2[1] - ca2[1]) * q]);
      var opp = Math.abs(((lo.k - hi.k) % 360 + 540) % 360 - 180) < 30;   /* реверсивное: противоположные направления */
      o.hours.push({ h: hh, v: Math.round(u2.v * 10) / 10, k: u2.v < 1e-6 ? lo.k : (opp ? (Math.cos((u2.k - lo.k) * R) > 0 ? lo.k : hi.k) : Math.round(u2.k)), tab: false, lo: lo, hi: hi });
    }
    var fl = num(d.cFl);
    if (!isFinite(fl)) { var mx = o.base.reduce(function (m, x) { return x.vc > m.vc ? x : m; }, o.base[0]); fl = mx.k; }
    o.fl = fl;
    o.hours.forEach(function (x) { x.s = (Math.cos((x.k - fl) * R) >= 0 ? 1 : -1) * x.v; });
    return o;
  }
  /* ---------- рисунки ---------- */
  function drawDepth(cv, LK, p) {
    var W = 1400, H = 880; cv.width = W; cv.height = H; var c = cv.getContext("2d"), FONT = LK.font, INK = LK.ink || "#1d3557";
    c.fillStyle = "#fff"; c.fillRect(0, 0, W, H);
    if (!p || !isFinite(p.Hz)) { c.fillStyle = "#555"; c.font = "24px " + FONT; c.fillText("Недостаточно данных для графика", 40, 60); return; }
    var ev = p.ev, t0 = Math.min(p.a.t, p.mvMin.t, p.T) - 1, t1 = Math.max(p.b.t, p.mvMin.t, p.T) + 1;
    var hs = ev.filter(function (e) { return e.t >= t0 - 8 && e.t <= t1 + 8; }).map(function (e) { return e.h; });
    var lev = p.hz + p.dhp, keel = lev - p.Ts;
    var top = Math.max.apply(null, hs.concat([lev])) + 1.3, bot = Math.min(-p.Hk, keel) - 0.8;
    var x0 = 540, x1 = W - 60, y0 = H - 70, y1 = 70;
    function X(t) { return x0 + (x1 - x0) * (t - t0) / (t1 - t0); } function Y(h) { return y0 - (y0 - y1) * (h - bot) / (top - bot); }
    /* грунт и нуль глубин */
    c.fillStyle = "#efe6d2"; c.fillRect(30, Y(-p.Hk), W - 60, Y(bot) - Y(-p.Hk) + 30);
    c.strokeStyle = "#8a6d3b"; c.lineWidth = 3; c.beginPath(); c.moveTo(30, Y(-p.Hk)); c.lineTo(W - 30, Y(-p.Hk)); c.stroke();
    c.font = "20px " + FONT; c.fillStyle = "#6b5326"; c.textAlign = "right"; c.fillText("грунт (Hк от нуля глубин)", W - 40, Y(-p.Hk) + 28);
    c.strokeStyle = "#777"; c.lineWidth = 1.5; c.setLineDash([10, 8]); c.beginPath(); c.moveTo(30, Y(0)); c.lineTo(W - 30, Y(0)); c.stroke(); c.setLineDash([]);
    c.fillStyle = "#555"; c.textAlign = "right"; c.fillText("0 глубин", W - 40, Y(0) - 8);
    /* кривая прилива */
    c.strokeStyle = INK; c.lineWidth = 4; c.beginPath(); var on = false;
    for (var k = 0; k <= 400; k++) { var tt = t0 + (t1 - t0) * k / 400, hh = hAt(ev, tt); if (!isFinite(hh)) { on = false; continue; } if (on) c.lineTo(X(tt), Y(hh)); else { c.moveTo(X(tt), Y(hh)); on = true; } }
    c.stroke();
    ev.forEach(function (e) { if (e.t < t0 || e.t > t1) return; c.fillStyle = e.pv ? "#c1121f" : "#2a9d8f"; c.beginPath(); c.arc(X(e.t), Y(e.h), 7, 0, 2 * Math.PI); c.fill();
      c.fillStyle = "#111"; c.font = "600 18px " + FONT; c.textAlign = "center"; c.fillText((e.pv ? "ПВ " : "МВ ") + T(e.t).slice(0, 5) + "  " + f(e.h, 1), Math.max(x0 + 60, Math.min(x1 - 60, X(e.t))), Y(e.h) + (e.pv ? -16 : 30)); });
    /* момент T и уровень */
    c.strokeStyle = "#999"; c.setLineDash([6, 6]); c.lineWidth = 2; c.beginPath(); c.moveTo(X(p.T), y1); c.lineTo(X(p.T), Y(-p.Hk)); c.stroke();
    c.strokeStyle = "#0b4f8a"; c.beginPath(); c.moveTo(30, Y(lev)); c.lineTo(X(p.T), Y(lev)); c.stroke(); c.setLineDash([]);
    c.fillStyle = "#0b4f8a"; c.beginPath(); c.arc(X(p.T), Y(p.hz), 7, 0, 2 * Math.PI); c.fill();
    c.fillStyle = "#333"; c.font = "600 20px " + FONT; c.textAlign = "center"; c.fillText("T = " + T(p.T).slice(0, 5), X(p.T), y1 - 16);
    c.fillStyle = "#0b4f8a"; c.font = "18px " + FONT; c.textAlign = "left"; c.fillText("уровень на T = hзад + Δhр", 360, Y(lev) - 10);
    /* судно */
    var sx = 250, bwT = 170, bwB = 104;
    c.fillStyle = "#e9eef5"; c.strokeStyle = "#222"; c.lineWidth = 3; c.beginPath();
    c.moveTo(sx - bwT / 2, Y(lev) - 34); c.lineTo(sx + bwT / 2, Y(lev) - 34); c.lineTo(sx + bwB / 2, Y(keel)); c.lineTo(sx - bwB / 2, Y(keel)); c.closePath(); c.fill(); c.stroke();
    c.beginPath(); c.rect(sx - 22, Y(lev) - 74, 44, 40); c.stroke();
    c.beginPath(); c.moveTo(sx, Y(lev) - 74); c.lineTo(sx, Y(lev) - 100); c.moveTo(sx - 18, Y(lev) - 94); c.lineTo(sx + 18, Y(lev) - 94); c.stroke();
    /* размеры: hзад (над нулём) и Hк (под нулём) — на одной вертикали справа от судна */
    G.dimV(c, 400, Y(0), Y(p.hz), "hзад " + f(p.hz, 2), INK, FONT, false);
    G.dimV(c, 400, Y(0), Y(-p.Hk), "Hк " + f(p.Hk, 1), "#6b5326", FONT, false);
    G.dimV(c, 140, Y(lev), Y(keel), "Tср " + f(p.Ts, 1), "#222", FONT, true);
    G.dimV(c, sx, Y(keel), Y(-p.Hk), "Hз " + f(p.Hz, 2), p.Hz >= 0 ? "#2f7d4f" : "#b53a2d", FONT, false);
    /* низкая МВ */
    c.strokeStyle = "#2a9d8f"; c.setLineDash([4, 6]); c.lineWidth = 2; c.beginPath(); c.moveTo(Math.max(x0, X(p.mvMin.t) - 120), Y(p.mvMin.h)); c.lineTo(Math.min(x1, X(p.mvMin.t) + 120), Y(p.mvMin.h)); c.stroke(); c.setLineDash([]);
    c.fillStyle = "#2a9d8f"; c.font = "600 19px " + FONT; c.textAlign = "center";
    c.fillText("Hмин = " + f(p.Hmin, 2) + " м (низкая МВ " + T(p.mvMin.t).slice(0, 5) + ")", Math.min(Math.max(X(p.mvMin.t), x0 + 170), x1 - 170), Y(p.mvMin.h) + 58);
    c.fillStyle = "#333"; c.font = "18px " + FONT; c.textAlign = "left"; c.fillText("Δhр = " + fs(p.dhp, 2) + " м;  B = " + f(p.B, 1) + " м", 40, H - 18);
    c.textAlign = "right"; c.fillText("T, ч", x1, H - 18);
  }
  function drawCur(cv, LK, o) {
    var W = 1400, H = 700; cv.width = W; cv.height = H; var c = cv.getContext("2d"), FONT = LK.font, INK = LK.ink || "#1d3557";
    c.fillStyle = "#fff"; c.fillRect(0, 0, W, H);
    var hs = (o && o.hours) || []; if (!hs.length) return;
    var mx = Math.max.apply(null, hs.map(function (x) { return Math.abs(x.s); }).concat([10])) / 10, top = Math.ceil(mx + 0.2), x0 = 110, x1 = W - 100, ym = H / 2, sy = (H / 2 - 70) / top;
    var ta = hs[0].h, tb = hs[hs.length - 1].h; if (tb === ta) tb = ta + 1;
    function X(t) { return x0 + (x1 - x0) * (t - ta) / (tb - ta); } function Y(v) { return ym - v * sy; }
    c.strokeStyle = "#e3e3e3"; c.lineWidth = 1; c.font = "19px " + FONT; c.fillStyle = "#444";
    for (var v = -top; v <= top + 1e-9; v += top > 4 ? 1 : 0.5) { c.beginPath(); c.moveTo(x0, Y(v)); c.lineTo(x1, Y(v)); c.stroke(); c.textAlign = "right"; c.fillText(f(Math.abs(v), top > 4 ? 0 : 1), x0 - 10, Y(v) + 7); }
    for (var t = ta; t <= tb; t++) { c.beginPath(); c.moveTo(X(t), Y(top)); c.lineTo(X(t), Y(-top)); c.stroke(); c.textAlign = "center"; c.fillText(String(t), X(t), ym + 26); }
    c.strokeStyle = "#333"; c.lineWidth = 2.5; c.beginPath(); c.moveTo(x0, ym); c.lineTo(x1, ym); c.stroke(); c.beginPath(); c.moveTo(x0, Y(top)); c.lineTo(x0, Y(-top)); c.stroke();
    c.fillStyle = "#333"; c.font = "600 20px " + FONT; c.textAlign = "left";
    c.fillText("V, уз — к " + deg(o.fl) + "° (прилив)", x0 + 10, Y(top) - 12); c.fillText("к " + deg(o.fl + 180) + "° (отлив)", x0 + 10, Y(-top) + 30);
    c.textAlign = "right"; c.fillText("T, ч", x1, ym - 12);
    c.strokeStyle = INK; c.lineWidth = 4; c.beginPath();
    var pts = hs.map(function (x) { return [X(x.h), Y(x.s / 10)]; });
    pts.forEach(function (pt, i) { if (!i) { c.moveTo(pt[0], pt[1]); return; } var p0 = pts[i - 1], mxp = (p0[0] + pt[0]) / 2; c.bezierCurveTo(mxp, p0[1], mxp, pt[1], pt[0], pt[1]); });
    c.stroke();
    hs.forEach(function (x, i) { c.fillStyle = x.tab ? "#c1121f" : "#2a9d8f"; c.beginPath(); c.arc(pts[i][0], pts[i][1], 7, 0, 2 * Math.PI); c.fill();
      c.fillStyle = "#111"; c.font = "600 17px " + FONT; c.textAlign = "center"; c.fillText(f(x.v / 10, 1) + " · " + deg(x.k) + "°", pts[i][0], pts[i][1] + (x.s >= 0 ? -14 : 28)); });
    c.font = "17px " + FONT; c.fillStyle = "#555"; c.textAlign = "left"; c.fillText("● из таблицы (чётные часы)   ● интерполяция", x0, H - 16);
  }
  function drawPlan(cv, LK, s) {
    var W = 1100, H = 1100; cv.width = W; cv.height = H; var c = cv.getContext("2d"), FONT = LK.font;
    c.fillStyle = "#fff"; c.fillRect(0, 0, W, H);
    if (!s || !isFinite(s.V)) return;
    var cx = W / 2, cy = H / 2, RR = 420, mx = Math.max(s.vq, s.vw, s.V, 0.5), stp = mx > 4 ? 1 : mx > 2 ? 0.5 : mx > 1 ? 0.25 : 0.1, nC = Math.ceil(mx / stp) + 1, sc = RR / (nC * stp);
    c.strokeStyle = "#d5d5d5"; c.lineWidth = 1;
    for (var a = 0; a < 360; a += 10) { var e = G.P2(cx, cy, a, RR); c.beginPath(); c.moveTo(cx, cy); c.lineTo(e[0], e[1]); c.stroke(); }
    c.font = "16px " + FONT; c.fillStyle = "#777"; c.textAlign = "left";
    for (var k = 1; k <= nC; k++) { c.beginPath(); c.arc(cx, cy, k * stp * sc, 0, 2 * Math.PI); c.stroke(); c.fillText(f(k * stp, stp < 0.25 ? 1 : 2).replace(/,?0+$/, "").replace(/,$/, ""), cx + 4, cy - k * stp * sc - 4); }
    G.compass(c, cx, cy, RR, FONT, 5, 10);
    var A = G.P2(cx, cy, s.kq, s.vq * sc), B = [A[0] + s.vw * sc * Math.sin(s.kw2 * R), A[1] - s.vw * sc * Math.cos(s.kw2 * R)], S = [B[0], B[1]];
    c.setLineDash([8, 7]); c.strokeStyle = "#aaa"; c.lineWidth = 2; var Bw = G.P2(cx, cy, s.kw2, s.vw * sc);
    c.beginPath(); c.moveTo(Bw[0], Bw[1]); c.lineTo(S[0], S[1]); c.stroke(); c.setLineDash([]);
    G.arrow(c, cx, cy, A[0], A[1], "#1d4e89", 4); G.arrow(c, A[0], A[1], B[0], B[1], "#2a9d8f", 4); G.arrow(c, cx, cy, Bw[0], Bw[1], "rgba(42,157,143,0.45)", 2.5, 14); G.arrow(c, cx, cy, S[0], S[1], "#c1121f", 5);
    c.font = "600 21px " + FONT; c.textAlign = "left";
    [["Vкв: " + deg(s.kq) + "° — " + f(s.vq, 2) + " уз", "#1d4e89"], ["Vт(в): " + deg(s.kw2) + "° — " + f(s.vw, 2) + " уз", "#2a9d8f"], ["VΣ: " + deg(s.K) + "° — " + f(s.V, 2) + " уз", "#c1121f"],
      ["Масштаб: 1 круг = " + String(stp).replace(".", ",") + " уз", "#444"]].forEach(function (x, i) { c.fillStyle = x[1]; c.fillText(x[0], 24, 34 + i * 30); });
  }

  App.taskWork({
    disc: "gmos", discName: "Гидрометеорологическое обеспечение судовождения", id: "gmos-pr3", order: 3, key: "gmos-pr3", file: "ГМОС_ПР3_КР3",
    no: "ПР 3", short: "ПР №3 (КР №3)", eyebrow: "ГМОС · Практическая работа №3 (Контрольная работа №3)", titleKind: "Контрольная работа",
    title: "Глубина под килем, приливо-отливное и суммарное течение",
    titleTopic: "Расчёт наименьшей ожидаемой глубины под килем и элементов суммарного течения",
    desc: "Два пункта: поправки к основному пункту, высота прилива на момент (закон косинуса), Hз и Hмин, графики. Приливо-отливное течение по N и C с интерполяцией, ветровое и суммарное течение, планшет.",
    variants: 1, variantData: function () { return Object.assign({}, DEF3); },
    varText: function (d) { return "Вариант " + (d.vno || "—"); },
    hint: "Выписывай с выданных страниц Таблиц приливов/течений: ПВ и МВ основного пункта на дату (строка «чч.мм высота»; воду предыдущих суток — с «-» в начале, следующих — с «+»), поправки из части II. Для пункта, который сам основной, поправки оставь пустыми. По умолчанию — данные варианта 21 и пример методички по течениям.",
    sections: [
      { title: "Работа", fields: [["vno", "Вариант №"], ["hsel", "Поправки высот", { select: [["sz", "сизигийные (как в методичке)"], ["kv", "квадратурные"]] }],
        ["P0", "Δhр: давление без поправки", { num: true, unit: "гПа" }], ["kP", "Δhр на 1 гПа", { num: true, unit: "м" }]] },
      { title: "Пункт № 1", fields: [["nm1", "Пункт", { wide: true }], ["T1", "Время (чч.мм)"], ["dt1", "Дата"], ["Hk1", "Нк", { num: true, unit: "м" }], ["P1", "P", { num: true, unit: "гПа" }], ["Ts1", "Тср (осадка)", { num: true, unit: "м" }],
        ["no1", "№ пункта"], ["co1", "Координаты"], ["op1", "№ ОП, стр."], ["tp1", "Поправка времени ПВ (±ч.мм)"], ["tm1", "Поправка времени МВ"],
        ["sp1", "сз. ПВ", { num: true, unit: "м" }], ["kp1", "кв. ПВ", { num: true, unit: "м" }], ["km1", "кв. МВ", { num: true, unit: "м" }], ["sm1", "сз. МВ", { num: true, unit: "м" }],
        ["ev1", "ОП на дату: «чч.мм высота» по строкам", { area: true, rows: 5 }]] },
      { title: "Пункт № 2", fields: [["nm2", "Пункт", { wide: true }], ["T2", "Время (чч.мм)"], ["dt2", "Дата"], ["Hk2", "Нк", { num: true, unit: "м" }], ["P2", "P", { num: true, unit: "гПа" }], ["Ts2", "Тср (осадка)", { num: true, unit: "м" }],
        ["no2", "№ пункта"], ["co2", "Координаты"], ["op2", "№ ОП, стр."], ["tp2", "Поправка времени ПВ (±ч.мм)"], ["tm2", "Поправка времени МВ"],
        ["sp2", "сз. ПВ", { num: true, unit: "м" }], ["kp2", "кв. ПВ", { num: true, unit: "м" }], ["km2", "кв. МВ", { num: true, unit: "м" }], ["sm2", "сз. МВ", { num: true, unit: "м" }],
        ["ev2", "ОП на дату: «чч.мм высота» по строкам", { area: true, rows: 5 }]] },
      { title: "2. Приливо-отливное течение", note: "Строки таблицы течений: «N V K V K …» — пары по выписанным чётным часам (V — в десятых узла, как в таблицах). Две строки (соседние N) — интерполяция по N.", fields: [
        ["cP", "Пункт", { wide: true }], ["cD", "Дата"], ["cT1", "С (ч)", { num: true }], ["cT2", "До (ч)", { num: true }], ["cN", "N (астр. данные)", { num: true }], ["cC", "C", { num: true }],
        ["cUseC", "Умножать на C", { select: [["1", "да (строка Vт × C)"], ["0", "нет (в таблицах сказано не вводить)"]] }], ["cPg", "Стр. ОП"], ["cFl", "Направление прилива (+), °", { num: true }],
        ["cH", "Часы таблицы (через пробел)", { wide: true }], ["cRows", "Строки «N V K …»", { area: true, rows: 3 }]] },
      { title: "3. Суммарное течение", fields: [["phi", "Широта φ", { num: true, unit: "°" }], ["hem", "Полушарие", { select: [["N", "северное (N)"], ["S", "южное (S)"]] }],
        ["kq", "Квазипост. течение: К", { num: true, unit: "°" }], ["vq", "V", { num: true, unit: "уз" }], ["kw", "Ветер: направление (откуда)", { num: true, unit: "°" }], ["w", "скорость W", { num: true, unit: "м/с" }],
        ["kc", "Коэфф. в U₀ = k·W/√sin φ", { num: true }]] }],
    goal: "Цель: рассчитать наименьшую ожидаемую глубину под килем корабля, стоящего на якоре, по Таблицам приливов; рассчитать элементы приливо-отливного течения по Таблицам течений, ветрового и суммарного течения.",
    figs: {
      d1: { caption: "Глубина под килем — пункт № 1", draw: function (cv, LK, d, r) { drawDepth(cv, LK, r.p[0]); } },
      d2: { caption: "Глубина под килем — пункт № 2", draw: function (cv, LK, d, r) { drawDepth(cv, LK, r.p[1]); } },
      cur: { caption: "Графическая интерполяция приливо-отливного течения", draw: function (cv, LK, d, r) { drawCur(cv, LK, r.cur); } },
      plan: { caption: "Построение суммарного течения на планшете", draw: function (cv, LK, d, r) { drawPlan(cv, LK, r.sum); } } },
    solve: function (d) {
      var warn = [], st = [], ps = [point(d, 1, warn), point(d, 2, warn)], sz = d.hsel !== "kv";
      st.push({ no: 1, noLabel: "1.1. ", title: "Таблица исходных данных", lines: [],
        tables: [{ headers: ["№", "Пункт", "Время", "Дата", "Нк, м", "P, гПа", "Тср, м"], rows: ps.map(function (p) { return [String(p.i), p.nm || "—", T(p.T), p.dt || "—", f(p.Hk, 1), f(p.P, 0), f(p.Ts, 1)]; }), widths: [0.8, 4.5, 2, 2, 1.8, 1.8, 1.8] }] });
      st.push({ no: 2, noLabel: "1.2. ", title: "Поправки к основному пункту", lines: ["Поправки высот берутся " + (sz ? "сизигийные (сз. ПВ, сз. МВ) — как в методичке" : "квадратурные (кв. ПВ, кв. МВ)") + "; для основного пункта поправок нет."],
        tables: [{ headers: ["№ пункта", "Координаты", "№ ОП, стр.", "ПВ", "МВ", "сз. ПВ", "кв. ПВ", "кв. МВ", "сз. МВ"],
          rows: ps.map(function (p) { function h(x) { return String(x || "").trim() ? fs(num(x), 1) : "—"; } return [p.no || String(p.i), p.co || "—", p.op || "—", dTtxt(p.tp), dTtxt(p.tm), h(p.sp), h(p.kp), h(p.km), h(p.sm)]; }), widths: [1.6, 3.2, 2, 1.4, 1.4, 1.3, 1.3, 1.3, 1.3], size: 9 }] });
      st.push({ no: 3, noLabel: "1.3–1.4. ", title: "Данные основного пункта и исправление поправками", lines: ps.map(function (p) { return "Пункт № " + p.i + ": T = TОП + ΔT (ПВ " + dTtxt(p.tp) + ", МВ " + dTtxt(p.tm) + ");  h = hОП + Δh (ПВ " + fs(p.cPV, 1) + ", МВ " + fs(p.cMV, 1) + ")"; }),
        tables: [].concat.apply([], ps.map(function (p) { return [
          { title: "1.3. Пункт № " + p.i + " — основной пункт" + (p.op ? " (" + p.op + ")" : ""), headers: ["ПВ: T", "h", "МВ: T", "h"], rows: evRows(p.ev0), widths: [3, 2, 3, 2] },
          { title: "1.4. Пункт № " + p.i + " — исправлено поправками", headers: ["ПВ: T", "h", "МВ: T", "h"], rows: evRows(p.ev), widths: [3, 2, 3, 2] }]; })) });
      function c5(fn) { return ps.map(function (p) { return isFinite(p.Hz) ? fn(p) : "—"; }); }
      st.push({ no: 4, noLabel: "1.5. ", title: "Расчёт времени падения/роста, величины прилива, глубины под килем", lines: [],
        after: [].concat.apply([], ps.filter(function (p) { return isFinite(p.Hz); }).map(function (p) { return [
          "Пункт № " + p.i + " (" + (p.nm || "") + ", " + T(p.T).slice(0, 5) + "): " + (p.rise ? "рост" : "падение") + " от " + (p.a.pv ? "ПВ " : "МВ ") + T(p.a.t) + " (" + f(p.a.h, 1) + ") до " + (p.b.pv ? "ПВ " : "МВ ") + T(p.b.t) + " (" + f(p.b.h, 1) + ")",
          "   T" + (p.rise ? "р" : "п") + " = " + T(p.b.t).slice(0, 5) + " − " + T(p.a.t).slice(0, 5) + " = " + dur(p.Tpr) + ";  B = " + f(p.pvW.h, 1) + " − " + f(p.mvW.h, 1) + " = " + f(p.B, 1) + " м;  ΔTбл = |" + T(p.T).slice(0, 5) + " − " + T(p.near.t).slice(0, 5) + "| = " + dur(p.dTb) + " (от " + (p.near.pv ? "ПВ" : "МВ") + ")",
          "   Δh = B/2·(1 − cos(180°·ΔTбл/T" + (p.rise ? "р" : "п") + ")) = " + f(p.B / 2, 2) + "·(1 − cos " + f(180 * p.dTb / p.Tpr, 1) + "°) = " + f(p.dh, 2) + " м (вспом. табл. 1);  hзад = " + f(p.near.h, 1) + (p.near.pv ? " − " : " + ") + f(p.dh, 2) + " = " + f(p.hz, 2) + " м",
          "   Δhр = (" + f(p.P0, 0) + " − " + f(p.P, 0) + ")·" + f(p.kP, 2) + " = " + fs(p.dhp, 2) + " м;  Hз = " + f(p.Hk, 1) + " + " + f(p.hz, 2) + " " + (p.dhp < 0 ? "− " : "+ ") + f(Math.abs(p.dhp), 2) + " − " + f(p.Ts, 1) + " = " + f(p.Hz, 2) + " м",
          "   Hмин (низкая МВ " + T(p.mvMin.t).slice(0, 5) + ", h = " + f(p.mvMin.h, 1) + ") = " + f(p.Hk, 1) + " + " + f(p.mvMin.h, 1) + " " + (p.dhp < 0 ? "− " : "+ ") + f(Math.abs(p.dhp), 2) + " − " + f(p.Ts, 1) + " = " + f(p.Hmin, 2) + " м"]; })),
        tables: [{ headers: ["Предварительные расчёты", "Формулы", "№ 1", "№ 2"], rows: [
          ["Время падения/роста", "Tп/р = Tмв/пв − Tпв/мв"].concat(c5(function (p) { return dur(p.Tpr) + (p.rise ? " (рост)" : " (падение)"); })),
          ["Величина прилива", "B = hПВ − hМВ"].concat(c5(function (p) { return f(p.B, 1) + " м"; })),
          ["Интервал от ближайшей ПВ или МВ", "ΔTбл = |T − Tмв/пв|"].concat(c5(function (p) { return dur(p.dTb) + " (" + (p.near.pv ? "ПВ " : "МВ ") + T(p.near.t).slice(0, 5) + ")"; })),
          ["Δh (вспом. табл. 1), интерполяция", "Δh = B/2·(1 − cos(180°·ΔT/Tп/р))"].concat(c5(function (p) { return (p.near.pv ? "−" : "+") + f(p.dh, 2); })),
          ["Высота прилива на Tзад", "hзад = hбл МВ/ПВ ± Δh"].concat(c5(function (p) { return f(p.hz, 2) + " м"; })),
          ["Поправка за давление Δhр", "Δhр = (1013 − P)·0,01"].concat(c5(function (p) { return fs(p.dhp, 2) + " м"; })),
          ["Глубина под килем Hз", "Hз = Hк + hзад + Δhр − Tср"].concat(c5(function (p) { return f(p.Hz, 2) + " м"; })),
          ["Мин. глубина под килем H", "Hмин = Hк + hМВ + Δhр − Tср"].concat(c5(function (p) { return f(p.Hmin, 2) + " м"; }))], widths: [4.2, 5.2, 3.6, 3.6], size: 9 }],
        answer: ps.map(function (p) { return "№ " + p.i + ": Hз = " + f(p.Hz, 2) + " м, Hмин = " + f(p.Hmin, 2) + " м"; }).join("; ") });
      st.push({ no: 5, noLabel: "1.6. ", title: "График", lines: ["Hз, Hк, hзад, Tср, B, hМВ — на графиках для каждого пункта."], figs: ["d1", "d2"] });
      /* течение */
      var cu = current(d, warn), cl = [];
      if (cu.base) {
        cl.push("Данные: " + (d.cP || "—") + ", " + (d.cD || "—") + ", с " + f(cu.t1, 0) + " до " + f(cu.t2, 0) + " ч;  N = " + f(cu.N, 1) + ", C = " + f(cu.C, 2) + (d.cPg ? ";  стр. " + d.cPg : ""));
        if (cu.r2) cl.push("Интерполяция по N между " + f(cu.r1.N, 1) + " и " + f(cu.r2.N, 1) + ": X = X₁ + (X₂ − X₁)·(" + f(cu.N, 1) + " − " + f(cu.r1.N, 1) + ")/(" + f(cu.r2.N, 1) + " − " + f(cu.r1.N, 1) + ") = X₁ + (X₂ − X₁)·" + f(cu.w, 2) + " (противоположные направления — со знаком)");
        else cl.push("Выписана одна строка N = " + f(cu.r1.N, 1) + " (ближайшая к N = " + f(cu.N, 1) + ") — интерполяция по N не требуется.");
        cl.push(cu.useC ? "Vт × C: скорость умножается на C = " + f(cu.C, 2) + "." : "Поправка за C не вводится.");
        cl.push("Нечётные часы — линейная интерполяция между соседними чётными; смена направления на противоположное — переход через 0.");
        var hdr = ["N"].concat(cu.H.map(function (h) { return h + " ч: Vт / Кт"; }));
        var trows = cu.rows.map(function (r) { return [f(r.N, 1)].concat(cu.H.map(function (h, j) { var x = r.val[j]; return x ? f(x.v, 0) + " / " + deg(x.k) : "—"; })); });
        if (cu.r2) trows.push([f(cu.N, 1)].concat(cu.base.map(function (x) { return f(x.v, 1) + " / " + deg(x.k); })));
        if (cu.useC) trows.push(["Vт × C"].concat(cu.base.map(function (x) { return f(x.vc, 1) + " / " + deg(x.k); })));
        st.push({ no: 6, noLabel: "2. ", title: "Расчёт приливо-отливного течения", lines: cl, figs: ["cur"],
          tables: [{ title: "2.2. Выбор курса и скорости течения (Vт — в десятых долях узла)", headers: hdr, rows: trows, size: 9 },
            { title: "Течение по часам (" + f(cu.t1, 0) + "–" + f(cu.t2, 0) + " ч)", headers: ["Час", "Vт, 0,1 уз", "V, уз", "Кт, °", "Источник"],
              rows: cu.hours.map(function (x) { return [String(x.h), f(x.v, 1), f(x.v / 10, 2), String(deg(x.k)), x.tab ? "таблица" : "интерполяция"]; }), widths: [1.5, 2.5, 2, 2, 3] }],
          answer: cu.hours.length ? cu.hours.map(function (x) { return x.h + " ч — " + deg(x.k) + "°, " + f(x.v / 10, 1) + " уз"; }).join("; ") : "" });
      }
      /* ветровое и суммарное */
      var phi = Math.abs(num(d.phi)), W = num(d.w), kw = num(d.kw), kq = num(d.kq), vq = num(d.vq), kc = num(d.kc) || 0.026, sh = d.hem === "S";
      var vw = Math.round(kc * W / Math.sqrt(Math.sin(phi * R)) * 100) / 100, kw2 = norm(kw + (sh ? 135 : 225));
      var ca = vec(vq, kq), cb = vec(vw, kw2), su = unvec([ca[0] + cb[0], ca[1] + cb[1]]);
      var sum = { vq: vq, kq: kq, vw: vw, kw2: kw2, K: su.k, V: Math.round(su.v * 100) / 100 };
      if ([phi, W, kw, kq, vq].every(isFinite)) {
        st.push({ no: 7, noLabel: "3.1. ", title: "Расчёт ветрового течения", lines: [
          "Kт(в) = Kw " + (sh ? "+ 135°" : "+ 225°") + " (" + (sh ? "южное полушарие — отклонение на 45° влево" : "северное полушарие — отклонение на 45° вправо") + ") = " + f(kw, 0) + (sh ? " + 135" : " + 225") + (kw + (sh ? 135 : 225) >= 360 ? " − 360" : "") + " = " + deg(kw2) + "°",
          "Vт(в) = U₀ = " + f(kc, 3) + "·W/√sin φ = " + f(kc, 3) + "·" + f(W, 0) + "/√sin " + f(phi, 0) + "° = " + f(kc * W, 3) + "/" + f(Math.sqrt(Math.sin(phi * R)), 3) + " = " + f(vw, 2) + " уз"],
          answer: "Kт(в) = " + deg(kw2) + "°, Vт(в) = " + f(vw, 2) + " уз" });
        st.push({ no: 8, noLabel: "3.2–3.3. ", title: "Суммарное течение", figs: ["plan"],
          lines: ["Векторное сложение квазипостоянного и ветрового течений (на планшете — треугольник векторов):",
            "Vx = " + f(vq, 2) + "·sin " + f(kq, 0) + "° + " + f(vw, 2) + "·sin " + deg(kw2) + "° = " + f(ca[0] + cb[0], 3) + ";  Vy = " + f(vq, 2) + "·cos " + f(kq, 0) + "° + " + f(vw, 2) + "·cos " + deg(kw2) + "° = " + f(ca[1] + cb[1], 3),
            "VтΣ = √(Vx² + Vy²) = " + f(sum.V, 2) + " уз;  KΣ = " + deg(su.k) + "°"],
          tables: [{ headers: ["φ", "К кв. теч.", "V кв. теч.", "Курс и скорость ветра", "Kт(в)", "Vт(в)", "KΣ", "VтΣ"],
            rows: [[f(phi, 0) + "° " + (sh ? "S" : "N"), deg(kq) + "°", f(vq, 1) + " уз", f(kw, 0) + "° — " + f(W, 0) + " м/с", deg(kw2) + "°", f(vw, 2) + " уз", deg(su.k) + "°", f(sum.V, 2) + " уз"]], widths: [1.6, 1.8, 1.8, 3.4, 1.6, 1.8, 1.6, 1.8], size: 9.5 }],
          answer: "KΣ = " + deg(su.k) + "°, VтΣ = " + f(sum.V, 2) + " уз" });
      } else warn.push("Суммарное течение: заполните φ, квазипостоянное течение и ветер.");
      return { steps: st, warn: warn, p: ps, cur: cu, sum: sum };
    }
  });
})();

/* ======================================================================
   МАНЕВРИРОВАНИЕ И УПРАВЛЕНИЕ СУДНОМ (МиУС)
   Методичка: Сливаев Б. Г. «Теоретические и практические основы управления
   судном», МГУ им. адм. Г. И. Невельского, 2027 («Методичка УС лабы 2027»).
   ЛР №1 «Влияние ветра на судно»:
     задание 1 — сила и момент ветра (ф-лы 1.1–1.4, табл. 1.1, прил. 2), графики;
     задание 2 — истинный ветер по табл. 1.3: круг СМО, графический способ,
                 аналитическая проверка.
   Исходник: папка «маневрирование/исходники приложения» (miyus1.js).
   ====================================================================== */
var MIYUS = (function () {
  "use strict";
  var DISC = App.discipline("miyus", "Маневрирование и управление судном", "МиУС");
  var R = Math.PI / 180;
  function num(x) { var s = String(x == null ? "" : x).trim().replace(/[−–]/g, "-").replace(",", "."); return s === "" ? NaN : parseFloat(s); }
  function rnd(x, n) { var p = Math.pow(10, n || 0), s = x < 0 ? -1 : 1; return s * Math.round(Math.abs(x) * p + 1e-7) / p; }
  function f(x, n) { return App.f(x, n === undefined ? 1 : n).replace(/^-/, "−"); }
  function fp(x, n) { var t = f(x, n); return t.charAt(0) === "−" ? "(" + t + ")" : t; }
  function fsg(x, n) { if (!isFinite(x)) return "—"; var t = f(x, n); return t.charAt(0) === "−" || /^0(,0+)?$/.test(t) ? t : "+" + t; }
  function deg(x) { x = ((x % 360) + 360) % 360; var r = Math.round(x); return r === 0 || r === 360 ? 360 : r; }
  function norm(x) { return ((x % 360) + 360) % 360; }
  function p3(x) { var s = String(Math.round(x)); while (s.length < 3) s = "0" + s; return s + "°"; }
  function g0(x) { return f(x, 1).replace(/,0$/, ""); }          /* 15 → «15», 7,5 → «7,5» */
  function lines(txt) { return String(txt || "").split(/\n+/).map(function (s) { return s.trim(); }).filter(function (s) { return s && !/^#/.test(s); }); }
  function toks(s) { return s.replace(/[−–]/g, "-").split(/[\s;|\/\t]+/).filter(Boolean); }
  function isNum(t) { return /^[+\-]?\d+(?:[.,]\d+)?$/.test(t); }

  /* ---------- прил. 2: тактико-технические данные моделей судов ---------- */
  var SHIPS = [
    { no: 1, name: "ALFA", gr: "α", Lmax: 120.3, Lpp: 105, B: 21.0, d: 7.5, disp: 13720, DW: 12349, Ba: 2250, Aa: 840 },
    { no: 2, name: "BETA", gr: "β", Lmax: 182.9, Lpp: 175, B: 24.4, d: 9.8, disp: 26771, DW: 22541, Ba: 2883, Aa: 781 },
    { no: 3, name: "GAMMA", gr: "γ", Lmax: 194.0, Lpp: 170, B: 30.5, d: 9.2, disp: 32000, DW: 22042, Ba: 3256, Aa: 946 },
    { no: 4, name: "DELTA", gr: "δ", Lmax: 174.0, Lpp: 168, B: 31.1, d: 12.0, disp: 55500, DW: 48765, Ba: 3300, Aa: 655 },
    { no: 5, name: "EPSILON", gr: "ε", Lmax: 330.7, Lpp: 320, B: 47.0, d: 11.1, disp: 135000, DW: 111773, Ba: 8463, Aa: 1880 },
    { no: 6, name: "DZETA", gr: "ζ", Lmax: 330.7, Lpp: 318.8, B: 46.5, d: 12.0, disp: 120000, DW: 106536, Ba: 8732, Aa: 1860 },
    { no: 7, name: "ETA", gr: "η", Lmax: 124.7, Lpp: 110.0, B: 20.0, d: 6.5, disp: 12989, DW: 10238, Ba: 3300, Aa: 560 },
    { no: 8, name: "TETA", gr: "θ", Lmax: 169.0, Lpp: 158.0, B: 27.2, d: 9.5, disp: 28080, DW: 22541, Ba: 5555, Aa: 920 },
    { no: 9, name: "YOTA", gr: "ι", Lmax: 190.0, Lpp: 178.0, B: 30.0, d: 10.0, disp: 33910, DW: 27315, Ba: 6940, Aa: 1000 },
    { no: 10, name: "KAPPA", gr: "κ", Lmax: 174.0, Lpp: 163.0, B: 29.8, d: 7.8, disp: 29990, DW: 25678, Ba: 4700, Aa: 841 },
    { no: 11, name: "LAMBDA", gr: "λ", Lmax: 289.0, Lpp: 279.0, B: 28.3, d: 13.4, disp: 86756, DW: 56732, Ba: 13950, Aa: 1350 },
    { no: 12, name: "MUE", gr: "μ", Lmax: 132.0, Lpp: 126.0, B: 20.0, d: 10.0, disp: 20241, DW: 16345, Ba: 3780, Aa: 600 }
  ];
  /* ---------- табл. 1.1: коэффициент силы ветра Ca (в грузу); столбец k — суда k и k+6 ---------- */
  var CA_PHI = [0, 10, 20, 30, 40, 50, 60, 90, 120, 150, 160, 170, 180];
  var CA = [
    [0.937, 1.169, 1.536, 1.712, 1.796, 1.482, 1.452, 1.343, 1.262, 1.516, 1.504, 0.946, 0.791],
    [0.684, 1.017, 1.374, 1.596, 1.599, 1.417, 1.355, 1.267, 1.306, 1.518, 1.494, 1.089, 0.721],
    [0.810, 0.919, 1.194, 1.364, 1.394, 1.414, 1.317, 1.176, 1.146, 1.402, 1.365, 1.061, 0.760],
    [0.778, 0.838, 1.171, 1.331, 1.365, 1.310, 1.200, 1.097, 1.197, 1.489, 1.417, 1.049, 0.705],
    [1.060, 1.024, 1.022, 1.016, 0.948, 0.973, 0.941, 0.963, 1.019, 1.261, 1.262, 1.133, 0.984],
    [0.815, 0.833, 0.984, 1.151, 1.228, 1.134, 1.064, 0.998, 0.961, 1.126, 1.113, 0.733, 0.578]
  ];
  /* ---------- табл. 1.3: элементы движения судна и кажущийся ветер по экипажам ---------- */
  var CREW = [
    [5, 10, "kn", 50, 10], [15, 7, "ms", 80, 4], [150, 9, "ms", 180, 10], [285, 7, "ms", 340, 11], [53, 12, "kn", 30, 12], [20, 6, "ms", 50, 10],
    [130, 8, "ms", 170, 9], [60, 14, "kn", 120, 9], [202, 7, "ms", 170, 8], [100, 6, "ms", 60, 5], [14, 12, "kn", 85, 8], [145, 7, "ms", 180, 13]
  ];
  function shipOf(v) { return SHIPS[(parseInt(v, 10) || 7) - 1] || SHIPS[6]; }
  function caText(v) { var col = ((parseInt(v, 10) || 7) - 1) % 6; return CA_PHI.map(function (p, i) { return p + " " + f(CA[col][i], 3); }).join("\n"); }
  function vd(v) {
    var s = shipOf(v), c = CREW[s.no - 1];
    return {
      variant: String(s.no), ship: s.name + " (" + s.gr + " " + s.no + ")",
      Lmax: f(s.Lmax, 1), Lpp: f(s.Lpp, 1), B: f(s.B, 1), d: f(s.d, 1), disp: String(s.disp), Ba: String(s.Ba), Aa: String(s.Aa),
      W1: "10", W2: "15", phis: "0 15 30 60 90 120 150 165 180", ca: caText(s.no), rho: "m", alf: "m", phiS: "30",
      K: String(c[0]), Vc: String(c[1]), vu: c[2], A: String(c[3]), Wa: String(c[4]), kn: "0,514"
    };
  }

  /* ---------- общие расчёты ---------- */
  function ptsOf(txt) { var p = []; lines(txt).forEach(function (l) { var t = toks(l).filter(isNum); if (t.length >= 2) p.push([num(t[0]), num(t[1])]); }); p.sort(function (a, b) { return a[0] - b[0]; }); return p; }
  function caAt(tab, x) {
    for (var i = 0; i < tab.length; i++) if (Math.abs(tab[i][0] - x) < 1e-9) return { v: tab[i][1], exact: true };
    if (tab.length < 2) return { v: tab.length ? tab[0][1] : NaN, exact: false, out: true };
    var k = 1; while (k < tab.length - 1 && x > tab[k][0]) k++;
    var a = tab[k - 1], b = tab[k];
    return { v: a[1] + (b[1] - a[1]) * (x - a[0]) / (b[0] - a[0]), a: a, b: b, exact: false, out: x < tab[0][0] || x > tab[tab.length - 1][0] };
  }
  function alphaOf(phi, corr) { var t = phi / 90; return (1 - 0.15 * (corr ? 1 - t : t) - 0.80 * Math.pow(1 - t, 3)) * 90; }
  function alphaTxt(phi, corr) {
    var p = g0(phi);
    return "[1 − 0,15·" + (corr ? "(1 − " + p + "/90)" : "(" + p + "/90)") + " − 0,80·(1 − " + p + "/90)³]·90";
  }
  function trueWind(K, V, A, W) {
    var x = W * Math.sin(A * R) - V * Math.sin(K * R), y = W * Math.cos(A * R) - V * Math.cos(K * R);
    var U = Math.sqrt(x * x + y * y), Ku = norm(Math.atan2(x, y) / R);
    var q = norm(A - K); if (q > 180) q -= 360;
    var Ua = Math.sqrt(V * V - 2 * W * V * Math.cos(q * R) + W * W);
    var cqi = Ua > 1e-9 ? (W * Math.cos(q * R) - V) / Ua : 1, qi = Math.acos(Math.max(-1, Math.min(1, cqi))) / R;
    return { x: x, y: y, U: U, K: Ku, q: q, Ua: Ua, cqi: cqi, qi: qi, sg: Math.sin(q * R) < 0 ? -1 : 1 };
  }
  function side(q) { return Math.abs(q) < 0.5 ? "прямо по носу" : Math.abs(Math.abs(q) - 180) < 0.5 ? "прямо по корме" : (q > 0 ? "правого" : "левого") + " борта"; }

  /* ---------- задание 1 ---------- */
  function task1(d, warn) {
    var Lpp = num(d.Lpp), Aa = num(d.Aa), Ba = num(d.Ba), W1 = num(d.W1), W2 = num(d.W2);
    var tech = d.rho === "t", rho = tech ? 0.000125 : 0.00125, corr = d.alf === "c", tab = ptsOf(d.ca);
    var phis = toks(String(d.phis || "")).filter(isNum).map(num).filter(function (p, i, a) { return a.indexOf(p) === i; }).sort(function (a, b) { return a - b; });
    if (![Lpp, Aa, Ba].every(isFinite)) warn.push("Задание 1: заполните Lpp, Aa и Ba (прил. 2).");
    if (![W1, W2].every(isFinite)) warn.push("Задание 1: заполните скорости ветра Wa1 и Wa2.");
    if (!tab.length) warn.push("Задание 1: не задана таблица Ca («φ Ca» по строкам, табл. 1.1).");
    if (!phis.length) warn.push("Задание 1: введите курсовые углы φ через пробел.");
    var rows = [];
    if (![Lpp, Aa, Ba, W1, W2].every(isFinite) || !tab.length) return { rows: rows, Lpp: Lpp, Aa: Aa, Ba: Ba, W1: W1, W2: W2, rho: rho, tech: tech, corr: corr };
    var half = rnd(Lpp / 2, 2);
    phis.forEach(function (p) {
      if (p < 0 || p > 180) { warn.push("φ = " + g0(p) + "° вне диапазона 0…180° — пропущен."); return; }
      var ci = caAt(tab, p), Ca = rnd(ci.v, 3);
      if (ci.out) warn.push("φ = " + g0(p) + "°: вне таблицы Ca — значение экстраполировано.");
      var c2 = rnd(Math.pow(Math.cos(p * R), 2), 4), s2 = rnd(Math.pow(Math.sin(p * R), 2), 4);
      var S = rnd(Aa * c2 + Ba * s2, 1);
      var Ra1 = rnd(0.5 * rho * Ca * S * W1 * W1, 2), Ra2 = rnd(0.5 * rho * Ca * S * W2 * W2, 2);
      var a = rnd((0.292 + 0.0023 * p) * Lpp, 2), al = rnd(alphaOf(p, corr), 1), CG = rnd(half - a, 2), sa = Math.sin(al * R);
      var Rm1 = rnd(Ra1 * sa * CG, 1), Rm2 = rnd(Ra2 * sa * CG, 1);
      rows.push({ p: p, ci: ci, Ca: Ca, c2: c2, s2: s2, S: S, Ra1: Ra1, Ra2: Ra2, a: a, al: al, CG: CG, sa: sa, Rm1: Rm1, Rm2: Rm2,
        Xa1: rnd(Ra1 * Math.cos(al * R), 2), Ya1: rnd(Ra1 * sa, 2) });
    });
    return { rows: rows, Lpp: Lpp, half: half, Aa: Aa, Ba: Ba, W1: W1, W2: W2, rho: rho, tech: tech, corr: corr };
  }
  function rhoTxt(t) { return t.tech ? "0,000125 тс·с²/м⁴" : "0,00125 т/м³"; }
  function t1Lines(t) {
    var L = ["Ra = ½·ρa·Ca·(Aa·cos²φ + Ba·sin²φ)·Wa²   (1.1),   ρa = " + rhoTxt(t) + ";  Aa = " + f(t.Aa, 0) + " м², Ba = " + f(t.Ba, 0) + " м²",
      "a = (0,292 + 0,0023·φ)·Lpp   (1.2),   Lpp = " + f(t.Lpp, 1) + " м",
      "α = [1 − 0,15·" + (t.corr ? "(1 − φ/90)" : "(φ/90)") + " − 0,80·(1 − φ/90)³]·90   (1.3)",
      "Rm = CG·Ya = Ra·sin α·(½Lpp − a)   (1.4),   ½Lpp = " + f(t.half, 2) + " м"];
    t.rows.forEach(function (r) {
      var p = g0(r.p), caS = r.ci.exact ? "Ca = " + f(r.Ca, 3) + " (табл. 1.1)" :
        "Ca = " + f(r.ci.a[1], 3) + " + (" + f(r.ci.b[1], 3) + " − " + f(r.ci.a[1], 3) + ")·(" + p + " − " + g0(r.ci.a[0]) + ")/(" + g0(r.ci.b[0]) + " − " + g0(r.ci.a[0]) + ") = " + f(r.Ca, 3) + " (интерполяция)";
      var k = "0,5·" + (t.tech ? "0,000125" : "0,00125") + "·" + f(r.Ca, 3) + "·" + f(r.S, 1);
      L.push("φ = " + p3(r.p) + ":  " + caS + ";  S = " + f(t.Aa, 0) + "·cos²" + p + "° + " + f(t.Ba, 0) + "·sin²" + p + "° = " + f(t.Aa, 0) + "·" + f(r.c2, 4) + " + " + f(t.Ba, 0) + "·" + f(r.s2, 4) + " = " + f(r.S, 1) + " м²");
      L.push("        Ra1 = " + k + "·" + g0(t.W1) + "² = " + f(r.Ra1, 2) + " тс;  Ra2 = " + k + "·" + g0(t.W2) + "² = " + f(r.Ra2, 2) + " тс");
      L.push("        a = (0,292 + 0,0023·" + p + ")·" + f(t.Lpp, 1) + " = " + f(r.a, 2) + " м;  α = " + alphaTxt(r.p, t.corr) + " = " + f(r.al, 1) + "°;  CG = " + f(t.half, 2) + " − " + f(r.a, 2) + " = " + f(r.CG, 2) + " м");
      L.push("        Rm1 = " + f(r.Ra1, 2) + "·sin " + f(r.al, 1) + "°·" + fp(r.CG, 2) + " = " + f(r.Rm1, 1) + " тс·м;  Rm2 = " + f(r.Ra2, 2) + "·sin " + f(r.al, 1) + "°·" + fp(r.CG, 2) + " = " + f(r.Rm2, 1) + " тс·м");
    });
    return L;
  }
  function t12(t) {
    var hd = ["φ"].concat(t.rows.map(function (r) { return p3(r.p); })), n = t.rows.length, w0 = 2.9, wi = n ? Math.min(1.8, (16.6 - w0) / n) : 1.5;
    function row(lab, k, nd) { return [lab].concat(t.rows.map(function (r) { return f(r[k], nd); })); }
    return { title: "Таблица 1.2 — Сила ветра и момент силы ветра", headers: hd,
      rows: [row("Ca", "Ca", 3), row("Ra1, тс (" + g0(t.W1) + " м/с)", "Ra1", 2), row("Ra2, тс (" + g0(t.W2) + " м/с)", "Ra2", 2), row("a, м", "a", 2), row("α, °", "al", 1), row("Rm1, тс·м", "Rm1", 1), row("Rm2, тс·м", "Rm2", 1)],
      widths: [w0].concat(t.rows.map(function () { return wi; })), size: n > 9 ? 8 : 9 };
  }
  function tAux(t) {
    return { title: "Вспомогательные величины", headers: ["φ", "cos²φ", "sin²φ", "S = Aa·cos²φ + Ba·sin²φ, м²", "½Lpp − a = CG, м", "sin α"],
      rows: t.rows.map(function (r) { return [p3(r.p), f(r.c2, 4), f(r.s2, 4), f(r.S, 1), f(r.CG, 2), f(r.sa, 4)]; }), widths: [1.6, 2, 2, 4.6, 3.6, 2], size: 9.5 };
  }
  function extremes(t, k) {
    var mx = null, mn = null;
    t.rows.forEach(function (r) { if (!mx || r[k] > mx[k]) mx = r; if (!mn || r[k] < mn[k]) mn = r; });
    return { mx: mx, mn: mn };
  }
  function t1Concl(t) {
    if (!t.rows.length) return [];
    var A = extremes(t, "Ra1"), M = extremes(t, "Rm1"), k = t.W2 / t.W1, out = [];
    var r0 = t.rows.filter(function (r) { return r.p === 0; })[0], r180 = t.rows.filter(function (r) { return r.p === 180; })[0];
    out.push("1. Сила ветра наибольшая при φ = " + p3(A.mx.p) + ": Ra1 = " + f(A.mx.Ra1, 2) + " тс, Ra2 = " + f(A.mx.Ra2, 2) + " тс — при ветре, близком к траверзному, работает проекция парусности на ДП (Ba = " + f(t.Ba, 0) + " м²), которая в " + f(t.Ba / t.Aa, 1) + " раза больше проекции на мидель (Aa = " + f(t.Aa, 0) + " м²). " +
      "Наименьшая сила — при ветре прямо по носу или по корме" + (r0 && r180 ? " (φ = 000°: " + f(r0.Ra1, 2) + " тс; φ = 180°: " + f(r180.Ra1, 2) + " тс при " + g0(t.W1) + " м/с)" : "") + ".");
    out.push("2. Сила и момент ветра пропорциональны квадрату скорости ветра: при увеличении скорости с " + g0(t.W1) + " до " + g0(t.W2) + " м/с (в " + f(k, 2).replace(/,?0+$/, "") + " раза) они возрастают в (" + g0(t.W2) + "/" + g0(t.W1) + ")² = " + f(k * k, 2).replace(/,?0+$/, "") + " раза.");
    var pos = t.rows.filter(function (r) { return r.Rm1 > 0; }), neg = t.rows.filter(function (r) { return r.Rm1 < 0; });
    if (pos.length) out.push("3. При носовых курсовых углах ветра (φ < 90°) точка приложения силы ветра C лежит в нос от центра тяжести (a < ½Lpp), момент Rm > 0 — ветер стремится увалить нос под ветер (развернуть судно от ветра). Наибольший момент — при φ = " + p3(M.mx.p) + ": Rm1 = " + f(M.mx.Rm1, 1) + " тс·м, Rm2 = " + f(M.mx.Rm2, 1) + " тс·м.");
    var r90 = t.rows.filter(function (r) { return r.p === 90; })[0];
    if (r90) out.push("4. При φ = 090° точка C почти совпадает с центром тяжести (a = " + f(r90.a, 2) + " м ≈ ½Lpp = " + f(t.half, 2) + " м), поэтому момент близок к нулю (Rm1 = " + f(r90.Rm1, 1) + " тс·м), хотя сила ветра почти максимальна: судно дрейфует, почти не разворачиваясь.");
    if (neg.length) out.push((r90 ? "5" : "4") + ". При кормовых курсовых углах (φ > 90°) точка C смещается в корму от центра тяжести, момент меняет знак (Rm < 0) — корма уваливается под ветер, нос приводится к ветру. Наибольший по модулю момент — при φ = " + p3(M.mn.p) + ": Rm1 = " + f(M.mn.Rm1, 1) + " тс·м, Rm2 = " + f(M.mn.Rm2, 1) + " тс·м.");
    return out;
  }

  /* ---------- задание 2 ---------- */
  function task2(d, warn) {
    var K = num(d.K), Vc = num(d.Vc), A = num(d.A), W = num(d.Wa), kn = num(d.kn) || 0.514, knots = d.vu === "kn";
    if (![K, Vc, A, W].every(isFinite)) { warn.push("Задание 2: заполните ИК, скорость судна, направление и скорость кажущегося ветра."); return null; }
    var V = knots ? rnd(Vc * kn, 1) : Vc, tw = trueWind(K, V, A, W);
    return { K: K, Vc: Vc, knots: knots, kn: kn, V: V, A: A, W: W, tw: tw, Ki: deg(tw.K), Vi: rnd(tw.U, 1) };
  }
  function t2Lines(x, crew) {
    var w = x.tw, L = [];
    L.push("Экипаж " + crew + ":  ИК = " + g0(x.K) + "°,  Vc = " + (x.knots ? g0(x.Vc) + " уз = " + g0(x.Vc) + "·" + f(x.kn, 3) + " = " + f(x.V, 1) + " м/с" : f(x.V, 1) + " м/с") +
      ";  кажущийся ветер: направление " + g0(x.A) + "°, Vв = " + g0(x.W) + " м/с");
    L.push("Аналитическая проверка:  q = " + g0(x.A) + " − " + g0(x.K) + " = " + fsg(w.q, 1).replace(/,0$/, "") + "° (КУ кажущегося ветра " + side(w.q) + ")");
    L.push("        Vи = √(Vc² + Vв² − 2·Vc·Vв·cos q) = √(" + f(x.V, 1) + "² + " + g0(x.W) + "² − 2·" + f(x.V, 1) + "·" + g0(x.W) + "·cos " + g0(Math.abs(w.q)) + "°) = " + f(w.Ua, 2) + " м/с");
    L.push("        cos qи = (Vв·cos q − Vc)/Vи = (" + g0(x.W) + "·cos " + g0(Math.abs(w.q)) + "° − " + f(x.V, 1) + ")/" + f(w.Ua, 2) + " = " + f(w.cqi, 4) + "  →  qи = " + f(w.qi, 1) + "° " + side(w.sg * w.qi));
    L.push("        Kи = ИК " + (w.sg < 0 ? "− " : "+ ") + "qи = " + g0(x.K) + (w.sg < 0 ? " − " : " + ") + f(w.qi, 1) + " = " + x.Ki + "°");
    return L;
  }
  function smoSteps(x) {
    var w = x.tw;
    return ["Масштаб: 1 деление сетки (1 см) = 1 м/с.",
      "1) Поворачиваем верхний круг так, чтобы деление " + g0(x.A) + "° (направление кажущегося ветра) встало у индекса.",
      "2) От центра по радиусу к индексу откладываем Vв = " + g0(x.W) + " м/с (" + g0(x.W) + " делений) — точка «В».",
      "3) Поворачиваем круг так, чтобы у индекса (указателя) встало деление ИК = " + g0(x.K) + "°.",
      "4) От центра по радиусу к указателю откладываем Vc = " + f(x.V, 1) + " м/с — точка «К».",
      "5) Вращаем круг, пока точки В и К не окажутся одна под другой на линии, параллельной диаметру через указатель (точка В ниже точки К).",
      "6) Против указателя отсчитываем направление истинного ветра Kи = " + x.Ki + "°; расстояние ВК = " + f(x.Vi, 1) + " делений — скорость истинного ветра Vи = " + f(x.Vi, 1) + " м/с."];
  }
  function grafSteps(x) {
    return ["1) Строим оси N–S и E–W.",
      "2) Из центра O по ИК = " + g0(x.K) + "° откладываем ОК = Vc = " + f(x.V, 1) + " м/с (в выбранном масштабе).",
      "3) Из центра O в том же масштабе по направлению кажущегося ветра " + g0(x.A) + "° откладываем ОВ = Vв = " + g0(x.W) + " м/с.",
      "4) Отрезок КВ, соединяющий концы векторов, переносим в центр — вектор ОИ: его направление Kи = " + x.Ki + "°, длина Vи = " + f(x.Vi, 1) + " м/с."];
  }

  /* ---------- контрольные вопросы ---------- */
  function questions(t, x) {
    var r90 = t.rows.filter(function (r) { return r.p === 90; })[0] || extremes(t, "Ra1").mx;
    var k2 = t.W2 && t.W1 ? f(Math.pow(t.W2 / t.W1, 2), 2).replace(/,?0+$/, "") : "2,25";
    return [
      { q: "Дать определение скорости ветра.",
        a: ["Скоростью ветра считается её среднее значение, измеренное за период более 10 минут на высоте 10 м. Порывы ветра могут превышать это среднее значение примерно в полтора раза; с увеличением высоты точки измерения средняя скорость ветра возрастает.",
          "Скорость ветра измеряют в м/с (узлах) или оценивают в баллах по 12-балльной шкале Бофорта, принятой ВМО; для ураганов используют шкалу Саффира — Симпсона (прил. 3)."] },
      { q: "Пояснить зависимость силы ветра от его скорости.",
        a: ["По уравнению Хьюза Ra = ½·ρa·Ca·(Aa·cos²φ + Ba·sin²φ)·Wa² сила ветра прямо пропорциональна квадрату его скорости: при увеличении скорости в n раз сила возрастает в n² раз. " +
          (r90 ? "В работе при увеличении скорости с " + g0(t.W1) + " до " + g0(t.W2) + " м/с сила при φ = " + p3(r90.p) + " выросла с " + f(r90.Ra1, 2) + " до " + f(r90.Ra2, 2) + " тс, т. е. в " + k2 + " раза; во столько же раз вырос и момент." : ""),
          "Кроме скорости, сила ветра зависит от площади парусности (проекций Aa и Ba), курсового угла кажущегося ветра (через cos²φ, sin²φ и коэффициент Ca, учитывающий обводы надводной части) и плотности воздуха."] },
      { q: "Вычислить направление и скорость истинного ветра, используя круг СМО и таблицу.",
        a: x ? ["Исходные данные (экипаж " + x.crew + "): ИК = " + g0(x.K) + "°, Vc = " + f(x.V, 1) + " м/с, кажущийся ветер " + g0(x.A) + "°, " + g0(x.W) + " м/с.",
          "На круге СМО: у индекса ставим " + g0(x.A) + "° и откладываем к индексу точку В (" + g0(x.W) + " м/с); ставим у указателя ИК " + g0(x.K) + "° и откладываем точку К (" + f(x.V, 1) + " м/с); вращаем круг, пока В и К не лягут на линию, параллельную диаметру через указатель (В ниже К). Против указателя — Kи = " + x.Ki + "°, длина ВК — Vи = " + f(x.Vi, 1) + " м/с.",
          "Расчёт (вместо таблицы) даёт то же: Vи = " + f(x.tw.Ua, 2) + " м/с, Kи = " + x.Ki + "° (см. задание 2)."] : ["См. задание 2."] },
      { q: "Объяснить влияние ветра на судно, лежащее в дрейфе.",
        a: ["Когда судно не имеет хода относительно воды, точка приложения силы ветра близка к середине корпуса, т. е. к центру тяжести судна, и плечо разворачивающего момента мало — ветер почти не разворачивает судно (в работе при φ = 090° момент близок к нулю).",
          "Обычно считают, что судно дрейфует по направлению ветра, если ветер действует по траверзу или несколько впереди или позади траверза. Скорость дрейфа тем больше, чем больше площадь парусности по сравнению с подводной частью корпуса."] },
      { q: "Как действует ветер траверзных направлений на судно, имеющее ход вперёд?",
        a: ["Точка приложения силы ветра остаётся вблизи середины корпуса, а центр вращения судна (Pivot Point) при ходе вперёд смещается в нос. Возникает плечо и вращающий момент, и нос судна приводится на ветер; величина момента зависит от силы ветра и скорости судна.",
          "Если снизить скорость, центр вращения смещается ещё дальше в нос (уменьшается гидродинамическая сила, препятствующая движению), плечо поворота и дрейф увеличиваются — тенденция усиливается. Это особенно важно при подходе к причалу с траверзным ветром, когда скорость приходится уменьшать. При ветре прямо по носу судно без дифферента управляется легко и не приводится на ветер вплоть до полной остановки."] }
    ];
  }

  