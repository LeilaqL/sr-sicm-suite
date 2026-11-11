/* Dual-bias triplet readout (dz99) - SR-SICM Suite v1.1
 * Leyla Nur Akay, Aarhus University, 2026.
 *
 * Port of sicm_tools/sicm_lib.py (read_park, load_curves, z_cross, blocks, triplets):
 *   z = Z detector x 0.879 nm/count, increasing away from the surface; channels 0 = Z command, 1 = Z detector, 4 = current
 *   I_inf = median current farther than `far` nm from the turnaround (charge default 1000 nm)
 *   z99 = crossing of I/I_inf = 0.99, from the bulk side (default) or from the surface side (option)
 *   blocks: consecutive curves with the same signed bias, < 5 min apart; the first 2 curves are discarded
 *   triplets: each block vs. the time-interpolated mean of its two opposite-polarity neighbours (same |V|)
 *   dz99 = z99(V-) - z99(V+); sign of dz99 = sign of the surface charge
 *   times from the TIFF header (tag 306, seconds); repetitions within a file are 14 s apart
 */
(function () {
  'use strict';

  const ZDET_NM = (typeof ZDET_DEFAULT !== 'undefined') ? ZDET_DEFAULT : 0.879;
  const NPTS = 4096, NCH = 5, CH_ZCMD = 0, CH_ZDET = 1, CH_CUR = 4, CUR_SCALE = 10 / 32768;
  const COL = { plus: '#d62728', minus: '#2f6fb5' };

  /* ------------------------------------------------------------ numeric core */
  const sortedCopy = (a) => Float64Array.from(a).sort();
  function pctSorted(s, p) {
    const n = s.length; if (!n) return NaN;
    const r = (p / 100) * (n - 1), i = Math.floor(r), f = r - i;
    return i + 1 < n ? s[i] + (s[i + 1] - s[i]) * f : s[i];
  }
  const median = (a) => pctSorted(sortedCopy(a), 50);
  const mean = (a) => a.reduce((s, x) => s + x, 0) / a.length;
  const sd1 = (a) => { if (a.length < 2) return NaN; const m = mean(a); return Math.sqrt(a.reduce((s, x) => s + (x - m) ** 2, 0) / (a.length - 1)); };

  /* crossing of level L. side 'bulk' (sicm_lib.z_cross): last point below L, interpolate to the next one.
   * side 'surface': first point at or above L counted from the surface, interpolate from the previous one. */
  function zCross(h, I, L, side) {
    const n = h.length;
    if (side === 'surface') {
      let j = -1;
      for (let i = 0; i < n; i++) if (I[i] >= L) { j = i; break; }
      if (j <= 0) return NaN;
      return h[j - 1] + (L - I[j - 1]) * (h[j] - h[j - 1]) / (I[j] - I[j - 1]);
    }
    let k = -1;
    for (let i = n - 1; i >= 0; i--) if (I[i] < L) { k = i; break; }
    if (k < 0 || k === n - 1) return NaN;
    return h[k] + (L - I[k]) * (h[k + 1] - h[k]) / (I[k + 1] - I[k]);
  }

  /* one approach curve from arrays in approach order (far -> turnaround), like sicm_lib.load_curves */
  function makeCurve(hApproach, Iabs, meta, opt) {
    const n = hApproach.length, hl = hApproach[n - 1];
    const far = [];
    for (let i = 0; i < n; i++) if (hApproach[i] - hl > opt.far) far.push(Iabs[i]);
    if (far.length < 20) return { skip: `curve shorter than far region (${opt.far} nm)` };
    const I0 = median(far);
    const idx = Array.from({ length: n }, (_, i) => i).sort((a, b) => hApproach[a] - hApproach[b] || a - b);
    const h = new Float64Array(n), I = new Float64Array(n);
    idx.forEach((j, i) => { h[i] = hApproach[j]; I[i] = Iabs[j] / I0; });
    const c = Object.assign({ h, I, I0, zturn: hl }, meta);
    c.z99 = zCross(h, I, opt.level, opt.side);
    return c;
  }

  function parseTiffTime(s) {
    const m = /^(\d{4}):(\d{2}):(\d{2}) (\d{2}):(\d{2}):(\d{2})/.exec(s || '');
    return m ? new Date(+m[1], +m[2] - 1, +m[3], +m[4], +m[5], +m[6]) : null;
  }

  /* sicm_lib.read_park + load_curves for one Park XE-Bio approach TIFF */
  function readParkCurves(buffer, name, opt) {
    const dv = new DataView(buffer), n = buffer.byteLength, skipped = [];
    const le = dv.getUint16(0) === 0x4949;
    if (!le && dv.getUint16(0) !== 0x4D4D) throw new Error('not a TIFF file');
    const ifd = dv.getUint32(4, le), cnt = dv.getUint16(ifd, le), tags = {};
    for (let i = 0; i < cnt; i++) {
      const o = ifd + 2 + 12 * i;
      tags[dv.getUint16(o, le)] = { count: dv.getUint32(o + 4, le), value: dv.getUint32(o + 8, le) };
    }
    if (!tags[50438] || !tags[50439]) throw new Error('no Park spectroscopy data (tags 50438/50439)');
    const zscale = Math.abs(dv.getFloat64(tags[50438].value + 80, le)) * 1000;     // nm per count, Z command
    let bias = NaN, setpoint = NaN, time = null;
    if (tags[50435]) { bias = dv.getFloat64(tags[50435].value + 212, le); setpoint = dv.getFloat64(tags[50435].value + 180, le); }
    if (tags[306] && tags[306].value + 19 <= n) {
      let s = ''; for (let i = 0; i < 19; i++) s += String.fromCharCode(dv.getUint8(tags[306].value + i));
      time = parseTiffTime(s);
    }
    const off = tags[50439].value, nval = Math.floor(tags[50439].count / 2);
    if (nval < NPTS * NCH) { skipped.push([name, 'short file (IV sweep?)']); return { curves: [], skipped }; }
    const nrep = Math.floor(nval / (NPTS * NCH));
    const at = (p, c, r) => dv.getInt16(off + 2 * (p + NPTS * (c + NCH * r)), le);   // numpy order 'F'
    const curves = [];
    for (let r = 0; r < nrep; r++) {
      let zmin = Infinity, zmax = -Infinity, kt = 0, sumI = 0;
      for (let p = 0; p < NPTS; p++) {
        const zc = at(p, CH_ZCMD, r) * zscale;
        if (zc > zmax) { zmax = zc; kt = p; }
        if (zc < zmin) zmin = zc;
        sumI += at(p, CH_CUR, r) * CUR_SCALE;
      }
      const label = nrep > 1 ? `${name} rep ${r + 1}` : name;
      if (zmax - zmin < 10) { skipped.push([label, 'empty repetition (Z command flat)']); continue; }
      const h = new Float64Array(kt + 1), I = new Float64Array(kt + 1);
      for (let p = 0; p <= kt; p++) { h[p] = -at(p, CH_ZDET, r) * opt.zdet; I[p] = Math.abs(at(p, CH_CUR, r) * CUR_SCALE); }
      const sgn = Math.sign(sumI);
      const V = isFinite(bias) && bias !== 0 ? bias : 0.1 * sgn;
      const t = time ? new Date(time.getTime() + r * opt.repS * 1000) : null;
      const c = makeCurve(h, I, { file: name, rep: r + 1, bias: V, setpoint, time: t }, opt);
      if (c.skip) { skipped.push([label, c.skip]); continue; }
      curves.push(c);
    }
    return { curves, skipped };
  }

  /* simulated demo (examples/triplet_demo/demo_curves.json) */
  function demoCurves(json, opt) {
    const curves = [], skipped = [];
    for (const d of json.curves) {
      const n = d.I_nA.length, h = new Float64Array(n), I = new Float64Array(n);
      for (let i = 0; i < n; i++) { h[i] = d.h0_nm + d.dh_nm * i; I[i] = Math.abs(d.I_nA[i]); }
      const t = new Date(d.time);
      const c = makeCurve(h, I, { file: d.file, rep: d.rep, bias: d.bias_V, setpoint: NaN, time: t }, opt);
      if (c.skip) skipped.push([d.file, c.skip]); else curves.push(c);
    }
    return { curves, skipped };
  }

  function sortCurves(C) {
    return C.slice().sort((a, b) => ((a.time ? a.time.getTime() : -Infinity) - (b.time ? b.time.getTime() : -Infinity)) || (a.file < b.file ? -1 : a.file > b.file ? 1 : 0));
  }

  /* sicm_lib.blocks */
  function blocks(curves, gapMin, discard) {
    const B = [];
    for (const c of curves) {
      const last = B.length ? B[B.length - 1] : null, lc = last ? last.curves[last.curves.length - 1] : null;
      if (last && Math.abs(last.bias - c.bias) < 1e-6 && c.time && lc.time && (c.time - lc.time) / 1000 < gapMin * 60) last.curves.push(c);
      else B.push({ bias: c.bias, curves: [c] });
    }
    for (const b of B) {
      b.used = b.curves.slice(discard);
      b.n = b.used.length;
      if (b.n) {
        b.t = mean(b.used.map((u) => u.time ? u.time.getTime() / 1000 : NaN)) / 60;
        const z = b.used.map((u) => u.z99).filter(isFinite);
        b.z99 = z.length ? mean(z) : NaN;                       // numpy nanmean
        const I0 = b.used.map((u) => u.I0);
        b.I0 = mean(I0);
        b.I0_drift_pct = 100 * (Math.max(...I0) - Math.min(...I0)) / b.I0;
      }
    }
    return B.filter((b) => b.n > 0);
  }

  /* sicm_lib.triplets */
  function triplets(B) {
    const out = [];
    for (let k = 1; k < B.length - 1; k++) {
      const a = B[k - 1], b = B[k], c = B[k + 1];
      if (Math.abs(a.bias - c.bias) > 1e-6 || Math.abs(a.bias - b.bias) < 1e-6) continue;
      const zi = a.z99 + (c.z99 - a.z99) * (b.t - a.t) / (c.t - a.t);
      const d = b.z99 - zi;
      out.push({ V: Math.round(Math.abs(a.bias - b.bias) / 2 * 1e4) / 1e4, dz99: b.bias < a.bias ? d : -d, t: b.t, k,
                 blocks: [a.curves[0].file, b.curves[0].file, c.curves[0].file],
                 maxDrift: Math.max(a.I0_drift_pct, b.I0_drift_pct, c.I0_drift_pct),
                 maxJump: Math.max(Math.abs(b.z99 - a.z99), Math.abs(c.z99 - b.z99)) });
    }
    return out;
  }

  function summarize(T, excludeFlagged) {
    const byV = new Map();
    for (const t of T) {
      if (excludeFlagged && t.maxDrift > 0.3) continue;
      if (!byV.has(t.V)) byV.set(t.V, []);
      byV.get(t.V).push(t.dz99);
    }
    return [...byV.entries()].sort((a, b) => a[0] - b[0]).map(([V, d]) => ({ V, n: d.length, mean: mean(d), sd: sd1(d) }));
  }

  /* ------------------------------------------------------------ UI */
  const S = { source: null, raw: null, demo: null, files: [], curves: [], skipped: [], B: [], T: [], expected: null };
  let root = null;
  const $ = (s) => root.querySelector(s);
  const esc = (s) => String(s).replace(/[&<>"']/g, (c) => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' }[c]));
  const f = (v, d = 2) => (isFinite(v) ? (v < 0 ? '−' : '') + Math.abs(v).toFixed(d) : '—');
  const fs = (v, d = 2) => (isFinite(v) ? (v > 0 ? '+' : v < 0 ? '−' : '') + Math.abs(v).toFixed(d) : '—');
  const SITE_ROOT = (() => { const s = document.currentScript && document.currentScript.src; return s ? s.replace(/js\/triplet\.js(\?.*)?$/, '') : ''; })();

  function opts() {
    const num = (id, d) => { const v = parseFloat($(id).value); return isFinite(v) ? v : d; };
    return { far: num('#trFar', 1000), level: num('#trLevel', 0.99), discard: Math.max(0, Math.round(num('#trDiscard', 2))),
             gapMin: num('#trGap', 5), repS: num('#trRep', 14), side: $('#trSide').value, zdet: ZDET_NM,
             excl: $('#trExcl').checked };
  }

  function build() {
    root.innerHTML = `
      <div class="cihd-head">
        <div>
          <h3 class="cihd-title">Dual-bias approach curves: triplet readout (Δz<sub>99</sub>)</h3>
          <p class="cihd-sub">Blocks of approach curves at alternating +V and −V. Each block is compared with the time-interpolated
            mean of its two opposite-polarity neighbours, which cancels linear drift. Δz<sub>99</sub> = z<sub>99</sub>(V−) − z<sub>99</sub>(V+);
            its sign is the sign of σ<sub>s</sub>. Converting Δz<sub>99</sub> to σ<sub>s</sub> needs a calibration for the same pipette, electrolyte and bias.</p>
        </div>
      </div>
      <div class="cihd-row">
        <button type="button" class="btn" id="trDemo">Load example (simulated demo)</button>
        <button type="button" class="btn" id="trOpen">Load Park TIFF files</button>
        <input type="file" id="trFiles" multiple accept=".tif,.tiff" hidden>
        <button type="button" class="btn" id="trClear">Clear</button>
      </div>
      <div class="cihd-permission" id="trDemoNote" hidden></div>
      <details class="cihd-step" open>
        <summary>Settings</summary>
        <div class="cihd-row">
          <label class="cihd-inline">I<sub>∞</sub> region (nm from turnaround) <input type="number" id="trFar" value="1000" step="50" min="50"></label>
          <label class="cihd-inline">Level <input type="number" id="trLevel" value="0.99" step="0.001"></label>
          <label class="cihd-inline">Discard first <input type="number" id="trDiscard" value="2" step="1" min="0"> curves per block</label>
        </div>
        <div class="cihd-row">
          <label class="cihd-inline">z<sub>99</sub> crossing
            <select id="trSide"><option value="bulk" selected>bulk side (default)</option><option value="surface">surface side (as in report)</option></select></label>
          <label class="cihd-inline">New block after <input type="number" id="trGap" value="5" step="1" min="1"> min</label>
          <label class="cihd-inline">Repetitions <input type="number" id="trRep" value="14" step="1" min="1"> s apart</label>
          <label class="cihd-check"><input type="checkbox" id="trExcl"> Exclude triplets with a flagged block (I∞ drift &gt; 0.3 %)</label>
        </div>
        <p class="cihd-note">Height: Z detector × ${ZDET_NM} nm/count. Times from the TIFF header (tag 306, seconds), not from the file name.</p>
      </details>
      <div class="cihd-msg" id="trMsg" role="status" aria-live="polite"></div>
      <div id="trOut" hidden>
        <div class="cihd-summary" id="trSummary"></div>
        <div class="tr-figs">
          <figure class="tr-fig"><div id="trFigBlocks"></div>
            <figcaption class="fig-cap">Block mean z<sub>99</sub> versus time (red +V, blue −V; x: min from the first block, y: nm). The drift between blocks is what the triplet scheme removes.</figcaption></figure>
          <figure class="tr-fig"><div id="trFigTrip"></div>
            <figcaption class="fig-cap">Δz<sub>99</sub> per triplet (nm). Open circles: triplets that use a flagged block. Dashed line: mean of the triplets included in the summary.</figcaption></figure>
        </div>
        <div class="cihd-note" id="trWarn"></div>
        <details class="cihd-step"><summary>Blocks</summary><div class="cihd-summary" id="trBlocks"></div></details>
        <details class="cihd-step"><summary>Triplets</summary><div class="cihd-summary" id="trTrips"></div></details>
        <details class="cihd-step"><summary>Skipped curves</summary><div class="cihd-note" id="trSkipped"></div></details>
        <div class="cihd-row cihd-export"><span class="cihd-label">Export</span>
          <button type="button" class="btn" id="trExpB">Blocks (CSV)</button><button type="button" class="btn" id="trExpT">Triplets (CSV)</button></div>
      </div>`;
    $('#trDemo').addEventListener('click', loadDemo);
    $('#trOpen').addEventListener('click', () => $('#trFiles').click());
    $('#trFiles').addEventListener('change', (e) => { loadFiles(e.target.files); e.target.value = ''; });
    $('#trClear').addEventListener('click', () => { Object.assign(S, { source: null, raw: null, demo: null, files: [], curves: [], B: [], T: [] }); $('#trOut').hidden = true; $('#trDemoNote').hidden = true; msg(''); });
    ['#trFar', '#trLevel', '#trDiscard', '#trSide', '#trGap', '#trRep'].forEach((id) => $(id).addEventListener('change', recompute));
    $('#trExcl').addEventListener('change', render);
    $('#trExpB').addEventListener('click', () => download('triplet_blocks.csv', blocksCSV()));
    $('#trExpT').addEventListener('click', () => download('triplet_triplets.csv', tripsCSV()));
    let rt = null;
    window.addEventListener('resize', () => { clearTimeout(rt); rt = setTimeout(() => { if (S.B.length) drawFigs(); }, 150); });
  }

  function msg(t, err) { const el = $('#trMsg'); el.textContent = t || ''; el.className = 'cihd-msg' + (err ? ' cihd-msg-err' : ''); }

  async function loadDemo() {
    try {
      msg('Loading simulated demo …');
      const [d, e] = await Promise.all([
        fetch(SITE_ROOT + 'examples/triplet_demo/demo_curves.json', { cache: 'no-cache' }).then((r) => { if (!r.ok) throw new Error('HTTP ' + r.status); return r.json(); }),
        fetch(SITE_ROOT + 'examples/triplet_demo/expected.json', { cache: 'no-cache' }).then((r) => r.ok ? r.json() : null).catch(() => null),
      ]);
      S.source = 'demo'; S.demo = d; S.expected = e; S.files = [];
      const note = $('#trDemoNote');
      note.hidden = false;
      note.innerHTML = `<strong>Simulated demo, not measured data.</strong> ${esc(d.description)}` +
        (e ? ` Known answer: Δz<sub>99</sub> = ${fs(e.true_dz99_nm, 1)} nm. Python reference (sicm_lib): ${fs(e.dz99_all.mean_nm)} ± ${f(e.dz99_all.sd_nm)} nm (n = ${e.dz99_all.n}).` : '');
      recompute();
    } catch (err) {
      msg((location.protocol === 'file:' ? 'The example needs a web server (http://, not file://). ' : 'Could not load the example: ') + err.message, true);
    }
  }

  async function loadFiles(list) {
    const files = Array.from(list || []);
    if (!files.length) return;
    msg(`Reading ${files.length} file(s) …`);
    const bufs = [];
    for (const fl of files) {
      try { bufs.push({ name: fl.name, buf: await fl.arrayBuffer() }); }
      catch (e) { bufs.push({ name: fl.name, error: e.message }); }
    }
    loadBuffers(bufs);
  }
  // [{name, buf}] -> analysis (also used by the tests)
  function loadBuffers(bufs) {
    S.source = 'files'; S.demo = null; S.expected = null; $('#trDemoNote').hidden = true;
    S.files = bufs.filter((b) => b.buf);
    recompute();
    for (const b of bufs.filter((x) => x.error)) S.skipped.push([b.name, b.error]);
  }

  function recompute() {
    if (!S.source) return;
    const o = opts();
    let C = [], sk = [];
    if (S.source === 'demo') ({ curves: C, skipped: sk } = demoCurves(S.demo, o));
    else for (const fl of S.files) {
      try { const r = readParkCurves(fl.buf, fl.name, o); C = C.concat(r.curves); sk = sk.concat(r.skipped); }
      catch (e) { sk.push([fl.name, e.message]); }
    }
    S.curves = sortCurves(C); S.skipped = sk;
    S.B = blocks(S.curves, o.gapMin, o.discard);
    S.T = triplets(S.B);
    msg(S.curves.length ? `${S.curves.length} curves, ${S.B.length} blocks, ${S.T.length} triplets.` : 'No usable approach curves.', !S.curves.length);
    render();
  }

  function render() {
    if (!S.source) return;
    const o = opts();
    $('#trOut').hidden = !S.B.length;
    if (!S.B.length) return;
    const all = summarize(S.T, false), sel = summarize(S.T, true);
    const rows = all.map((a) => {
      const b = sel.find((x) => x.V === a.V) || { n: 0, mean: NaN, sd: NaN };
      return `<tr><td>±${Math.round(a.V * 1000)} mV</td><td class="num">${a.n}</td><td class="num">${fs(a.mean)} ± ${f(a.sd)}</td>
              <td class="num">${b.n}</td><td class="num">${fs(b.mean)} ± ${f(b.sd)}</td></tr>`;
    }).join('');
    const pick = o.excl ? 'without flagged' : 'all triplets';
    $('#trSummary').innerHTML = `<table class="cihd-table"><caption>Δz<sub>99</sub> per |V| (mean ± SD, nm; z<sub>99</sub> from the ${o.side} side). Shown in the figure: ${pick}.</caption>
      <thead><tr><th>|V|</th><th>n</th><th>All triplets</th><th>n</th><th>Without flagged blocks</th></tr></thead><tbody>${rows ||
      '<tr><td colspan="5">No triplets: the series needs at least three blocks alternating +V / −V at the same |V|.</td></tr>'}</tbody></table>`;
    // warnings (sicm_lib conventions)
    const w = [];
    const fl = S.B.map((b, i) => [i + 1, b]).filter(([, b]) => b.I0_drift_pct > 0.3);
    if (fl.length) w.push(`I<sub>∞</sub> varies by more than 0.3 % within block ${fl.map(([i, b]) => `${i} (${b.I0_drift_pct.toFixed(2)} %)`).join(', ')}; results are given with and without the triplets that use it.`);
    const jumps = S.T.filter((t) => t.maxJump > 500);
    if (jumps.length) w.push(`${jumps.length} triplet(s) span a z<sub>99</sub> jump &gt; 500 nm between blocks (re-engagement?); such triplets are not valid.`);
    const nanB = S.B.filter((b) => !isFinite(b.z99)).length;
    if (nanB) w.push(`${nanB} block(s) without a level crossing.`);
    if (S.source === 'files' && S.curves.some((c) => !c.time)) w.push('Some files have no time stamp (tag 306); blocks cannot be formed for them.');
    $('#trWarn').innerHTML = w.length ? '<strong>Warnings.</strong> ' + w.join(' ') : '';
    // tables
    const td = (v) => `<td class="num">${v}</td>`;
    $('#trBlocks').innerHTML = `<table class="cihd-table"><thead><tr><th>#</th><th>First file</th><th>Bias (V)</th><th>n used</th><th>z<sub>99</sub> (nm)</th><th>I<sub>∞</sub> (nA)</th><th>I<sub>∞</sub> drift (%)</th></tr></thead><tbody>` +
      S.B.map((b, i) => `<tr><td>${i + 1}</td><td>${esc(b.curves[0].file)}</td>${td(fs(b.bias, 3))}${td(b.n)}${td(f(b.z99, 1))}${td(f(b.I0, 4))}${td(f(b.I0_drift_pct, 2) + (b.I0_drift_pct > 0.3 ? ' &gt; 0.3' : ''))}</tr>`).join('') + '</tbody></table>';
    $('#trTrips').innerHTML = `<table class="cihd-table"><thead><tr><th>Middle block</th><th>|V| (mV)</th><th>Δz<sub>99</sub> (nm)</th><th>max I<sub>∞</sub> drift (%)</th></tr></thead><tbody>` +
      S.T.map((t) => `<tr><td>${t.k + 1}: ${esc(t.blocks[1])}</td>${td(Math.round(t.V * 1000))}${td(fs(t.dz99))}${td(f(t.maxDrift, 2))}</tr>`).join('') + '</tbody></table>';
    $('#trSkipped').innerHTML = S.skipped.length ? S.skipped.map(([a, b]) => `${esc(a)}: ${esc(b)}`).join('<br>') : 'None.';
    drawFigs();
  }

  /* ------------------------------------------------------------ figures */
  const W = (el) => Math.max(300, Math.min(560, Math.round(el.clientWidth || 520)));
  function ticks(lo, hi, want) {
    const raw = (hi - lo) / Math.max(1, want), mag = Math.pow(10, Math.floor(Math.log10(raw || 1)));
    const step = [1, 2, 2.5, 5, 10].map((m) => m * mag).find((s) => s >= raw) || mag * 10;
    const out = []; for (let v = Math.ceil(lo / step) * step; v <= hi + 1e-9; v += step) out.push(+v.toFixed(10)); return out;
  }
  function axes(w, h, m, x0, x1, y0, y1, xl, yl) {
    const X = (v) => m.l + ((v - x0) / (x1 - x0)) * (w - m.l - m.r), Y = (v) => m.t + (1 - (v - y0) / (y1 - y0)) * (h - m.t - m.b);
    let s = `<svg viewBox="0 0 ${w} ${h}" xmlns="http://www.w3.org/2000/svg" font-family="Geist, Helvetica, Arial, sans-serif">`;
    for (const t of ticks(y0, y1, 5)) s += `<line x1="${m.l}" x2="${w - m.r}" y1="${Y(t)}" y2="${Y(t)}" stroke="var(--line-2)"/><text x="${m.l - 6}" y="${Y(t) + 4}" font-size="11" text-anchor="end" fill="var(--faint)">${String(t).replace('-', '−')}</text>`;
    for (const t of ticks(x0, x1, 6)) s += `<text x="${X(t)}" y="${h - m.b + 16}" font-size="11" text-anchor="middle" fill="var(--faint)">${t}</text>`;
    s += `<line x1="${m.l}" x2="${m.l}" y1="${m.t}" y2="${h - m.b}" stroke="var(--muted)"/><line x1="${m.l}" x2="${w - m.r}" y1="${h - m.b}" y2="${h - m.b}" stroke="var(--muted)"/>`;
    s += `<text x="${(m.l + w - m.r) / 2}" y="${h - 6}" font-size="11" text-anchor="middle" fill="var(--faint)">${xl}</text>`;
    s += `<text transform="translate(12 ${(m.t + h - m.b) / 2}) rotate(-90)" font-size="11" text-anchor="middle" fill="var(--faint)">${yl}</text>`;
    return { s, X, Y };
  }
  function drawFigs() {
    const o = opts(), m = { l: 58, r: 12, t: 12, b: 40 }, H = 240;
    // blocks
    const el1 = $('#trFigBlocks'), w1 = W(el1), B = S.B.filter((b) => isFinite(b.z99) && isFinite(b.t));
    if (B.length) {
      const t0 = B[0].t, xs = B.map((b) => b.t - t0), ys = B.map((b) => b.z99);
      const pad = 0.08 * ((Math.max(...ys) - Math.min(...ys)) || 10);
      const A = axes(w1, H, m, Math.min(0, ...xs), Math.max(...xs) * 1.04 + 0.1, Math.min(...ys) - pad, Math.max(...ys) + pad, 'time (min)', 'block z99 (nm)');
      let s = A.s;
      B.forEach((b, i) => {
        const c = b.bias > 0 ? COL.plus : COL.minus, fl = b.I0_drift_pct > 0.3;
        s += `<circle cx="${A.X(xs[i])}" cy="${A.Y(ys[i])}" r="5" fill="${fl ? 'var(--surface)' : c}" stroke="${c}" stroke-width="2"><title>Block ${S.B.indexOf(b) + 1}, ${b.bias > 0 ? '+' : '−'}${Math.abs(b.bias * 1000).toFixed(0)} mV: z99 ${b.z99.toFixed(1)} nm, I∞ drift ${b.I0_drift_pct.toFixed(2)} %</title></circle>`;
      });
      const lx = w1 - m.r - 78;   // legend top right (drift usually runs downwards to the right)
      s += `<g font-size="11" fill="var(--txt)"><circle cx="${lx}" cy="${m.t + 10}" r="4.5" fill="${COL.plus}"/><text x="${lx + 9}" y="${m.t + 14}">+V</text>` +
           `<circle cx="${lx + 40}" cy="${m.t + 10}" r="4.5" fill="${COL.minus}"/><text x="${lx + 49}" y="${m.t + 14}">−V</text></g>`;
      el1.innerHTML = s + '</svg>';
    } else el1.innerHTML = '';
    // triplets
    const el2 = $('#trFigTrip'), w2 = W(el2), T = S.T;
    if (T.length) {
      const ys = T.map((t) => t.dz99), inc = T.filter((t) => !(o.excl && t.maxDrift > 0.3)).map((t) => t.dz99);
      const lo = Math.min(0, ...ys), hi = Math.max(0, ...ys), pad = 0.12 * ((hi - lo) || 1);
      const A = axes(w2, H, m, 0.5, T.length + 0.5, lo - pad, hi + pad, 'triplet', 'Δz99 (nm)');
      let s = A.s + `<line x1="${m.l}" x2="${w2 - m.r}" y1="${A.Y(0)}" y2="${A.Y(0)}" stroke="var(--muted)"/>`;
      if (inc.length) { const mu = mean(inc); s += `<line x1="${m.l}" x2="${w2 - m.r}" y1="${A.Y(mu)}" y2="${A.Y(mu)}" stroke="var(--acc)" stroke-dasharray="6 4" stroke-width="1.5"/><line x1="${m.l + 8}" x2="${m.l + 26}" y1="${m.t + 10}" y2="${m.t + 10}" stroke="var(--acc)" stroke-dasharray="6 4" stroke-width="1.5"/><text x="${m.l + 31}" y="${m.t + 14}" font-size="11" fill="var(--txt)">mean ${fs(mu)} nm (n = ${inc.length})</text>`; }
      T.forEach((t, i) => {
        const fl = t.maxDrift > 0.3;
        s += `<circle cx="${A.X(i + 1)}" cy="${A.Y(t.dz99)}" r="5" fill="${fl ? 'var(--surface)' : 'var(--txt)'}" stroke="var(--txt)" stroke-width="1.6"><title>Triplet ${i + 1} (middle block ${t.k + 1}): Δz99 ${fs(t.dz99)} nm${fl ? ', uses a flagged block' : ''}</title></circle>`;
      });
      el2.innerHTML = s + '</svg>';
    } else el2.innerHTML = '';
  }

  /* ------------------------------------------------------------ export */
  function download(name, text) {
    const a = document.createElement('a'); a.href = URL.createObjectURL(new Blob([text], { type: 'text/csv' })); a.download = name;
    document.body.appendChild(a); a.click(); a.remove(); setTimeout(() => URL.revokeObjectURL(a.href), 2000);
  }
  function header() { const o = opts(); return `# SR-SICM Suite v1.1 triplet readout; ${S.source === 'demo' ? 'SIMULATED DEMO DATA; ' : ''}far ${o.far} nm, level ${o.level}, ${o.side} side, discard ${o.discard}, Z detector ${o.zdet} nm/count\n`; }
  function blocksCSV() { return header() + 'block,first_file,bias_V,n_used,z99_nm,I_inf_nA,I_inf_drift_pct\n' + S.B.map((b, i) => [i + 1, b.curves[0].file, b.bias, b.n, b.z99.toFixed(3), b.I0.toFixed(5), b.I0_drift_pct.toFixed(3)].join(',')).join('\n') + '\n'; }
  function tripsCSV() { return header() + 'V,dz99_nm,block_prev,block_mid,block_next,max_I_inf_drift_pct\n' + S.T.map((t) => [t.V, t.dz99.toFixed(3), ...t.blocks, t.maxDrift.toFixed(3)].join(',')).join('\n') + '\n'; }

  function init() { if (root) return true; root = document.getElementById('trip-root'); if (root) build(); return !!root; }
  if (document.readyState === 'loading') document.addEventListener('DOMContentLoaded', init); else init();
  const ready = (fn) => (...a) => (init() ? fn(...a) : undefined);   // public calls work before DOMContentLoaded too

  window.TRIPLET = { core: { zCross, makeCurve, readParkCurves, demoCurves, sortCurves, blocks, triplets, summarize, median },
                     loadDemo: ready(loadDemo), loadFiles: ready(loadFiles), loadBuffers: ready(loadBuffers), state: S,
                     set: ready((id, v) => { const el = $(id); if (el.type === 'checkbox') el.checked = v; else el.value = v; el.dispatchEvent(new Event('change')); }) };
})();
