/* SR-SICM Suite v1.1 - shared page content: conventions box and the two "what the outputs mean" figures.
 * Leyla Nur Akay, Aarhus University, 2026. */
(function () {
  'use strict';

  // One list of conventions, shown wherever a .conv-slot exists. Numbers are taken from the code where possible.
  const ZDET = (typeof ZDET_DEFAULT !== 'undefined') ? ZDET_DEFAULT : 0.879;
  const CONV = [
    ['Δz<sub>99</sub>', 'z<sub>99</sub>(V−) − z<sub>99</sub>(V+), from approach curves. The sign of Δz<sub>99</sub> is the sign of the surface charge σ<sub>s</sub>. ' +
      'Called Δh in Klausen et al. (2016).'],
    ['z<sub>99</sub>', 'Position where I/I<sub>∞</sub> crosses 0.99. Default: approached from the bulk side. The triplet and stiffness modules also offer ' +
      '"surface side" (first crossing counted from the surface) as an option.'],
    ['I<sub>∞</sub>', 'Median current farther than a set distance from the turnaround: 1000 nm for charge (long curves), 300 nm for stiffness.'],
    ['Height', `Z detector × ${ZDET} nm/count (calibration on a PDMS grating, 3 Nov 2026).`],
    ['Slope s', 'Linear fit of I/I<sub>∞</sub> versus z between 0.98 and 0.99, one fit per curve (never on an averaged curve).'],
    ['Discarded curves', 'Charge: the first 2 curves of each block. Stiffness: the first curve of each file.'],
    ['Blocks', 'Consecutive curves at one bias, normally 5. A block (charge) or file (stiffness) whose I<sub>∞</sub> varies by more than 0.3 % is flagged.'],
    ['Triplets', 'Each block is compared with the time-interpolated mean of its two opposite-polarity neighbours, which cancels linear drift. ' +
      'Times from the TIFF header (tag 306, with seconds); repetitions within a file are 14 s apart.'],
    ['Stiffness', 's<sub>∞</sub>/s and indentation δ = 0.01 (1/s − 1/s<sub>∞</sub>), with s<sub>∞</sub> from a rigid reference measured with the same ' +
      'pipette and the same bias. E/p<sub>0</sub> is an apparent value; absolute E (Pa) is not reported.'],
    ['Δh (imaging)', 'Δh = h(−V) − h(+V) from two topography images, on the overlap after registration. Δh (imaging) follows the same sign rule as Δz<sub>99</sub>.'],
  ];
  function fillConventions() {
    const html = '<div class="conv"><h3>Conventions used on this site</h3><dl>' +
      CONV.map(([k, v]) => `<dt>${k}</dt><dd>${v}</dd>`).join('') +
      '</dl><p class="conv-note">Colours: +V red, −V blue, Δh green.</p></div>';
    document.querySelectorAll('.conv-slot').forEach((el) => { el.innerHTML = html; });
  }

  const NS = 'http://www.w3.org/2000/svg';
  const fmt = (v, d = 1) => (v > 0 ? '+' : v < 0 ? '−' : '') + Math.abs(v).toFixed(d);

  // Box 1: median dh on lipid patches (CIHD example, published re-processed values)
  // drawing width follows the box (320-460 user units), so text keeps ~11 px on phones
  const boxW = (box) => Math.max(320, Math.min(460, Math.round(box.clientWidth || 460)));
  let LIPIDS = null;
  async function chargeFigure() {
    const box = document.getElementById('useChargeFig');
    if (!box) return;
    if (!LIPIDS) {
      try {
        const r = await fetch('examples/cihd/manifest.json', { cache: 'no-cache' });
        LIPIDS = (await r.json()).lipids;
      } catch (e) {
        LIPIDS = [{ id: 'DPTAP', charge: 'cationic', published_dh_nm: { median: 9.2, q25: 6.2, q75: 11.3 } },
                  { id: 'DPPE', charge: 'zwitterionic', published_dh_nm: { median: 7.0, q25: 6.5, q75: 7.4 } },
                  { id: 'DPPG', charge: 'anionic', published_dh_nm: { median: -2.1, q25: -3.4, q75: -1.2 } }];
      }
    }
    const L = LIPIDS;
    const W = boxW(box), H = 250, ml = 46, mr = 12, mt = 16, mb = 46, pw = W - ml - mr, ph = H - mt - mb;
    const lo = -6, hi = 14;
    const Y = (v) => mt + (1 - (v - lo) / (hi - lo)) * ph;
    const bw = pw / L.length;
    let s = `<svg viewBox="0 0 ${W} ${H}" xmlns="${NS}" font-family="Geist, Helvetica, Arial, sans-serif">`;
    for (const t of [-5, 0, 5, 10]) {
      s += `<line x1="${ml}" x2="${ml + pw}" y1="${Y(t)}" y2="${Y(t)}" stroke="var(--line-2)" stroke-width="1"/>`;
      s += `<text x="${ml - 6}" y="${Y(t) + 4}" font-size="11" text-anchor="end" fill="var(--faint)">${t === 0 ? '0' : fmt(t, 0)}</text>`;
    }
    s += `<line x1="${ml}" x2="${ml + pw}" y1="${Y(0)}" y2="${Y(0)}" stroke="var(--muted)" stroke-width="1"/>`;
    L.forEach((l, i) => {
      const p = l.published_dh_nm, cx = ml + bw * (i + 0.5), w = Math.min(46, bw * 0.42);
      const col = p.median >= 0 ? '#d62728' : '#2f6fb5';
      const y0 = Y(0), y1 = Y(p.median), top = Math.min(y0, y1), h = Math.max(1, Math.abs(y1 - y0));
      s += `<g><title>${l.id} (${l.charge}): median Δh ${fmt(p.median)} nm, IQR ${fmt(p.q25)} to ${fmt(p.q75)} nm</title>`;
      s += `<rect x="${cx - w / 2}" y="${top}" width="${w}" height="${h}" fill="${col}" rx="2"/>`;
      s += `<line x1="${cx}" x2="${cx}" y1="${Y(p.q25)}" y2="${Y(p.q75)}" stroke="var(--txt)" stroke-width="1.2"/>`;
      s += `<line x1="${cx - 6}" x2="${cx + 6}" y1="${Y(p.q25)}" y2="${Y(p.q25)}" stroke="var(--txt)" stroke-width="1.2"/>`;
      s += `<line x1="${cx - 6}" x2="${cx + 6}" y1="${Y(p.q75)}" y2="${Y(p.q75)}" stroke="var(--txt)" stroke-width="1.2"/>`;
      const ly = p.median >= 0 ? Math.min(Y(p.q75), top) - 6 : Math.max(Y(p.q25), top + h) + 14;
      s += `<text x="${cx}" y="${ly}" font-size="12" font-weight="600" text-anchor="middle" fill="var(--txt)">${fmt(p.median)} nm</text>`;
      s += `<text x="${cx}" y="${H - mb + 18}" font-size="12" text-anchor="middle" fill="var(--txt)">${l.id}</text>`;
      s += `<text x="${cx}" y="${H - mb + 32}" font-size="11" text-anchor="middle" fill="var(--faint)">${l.charge}</text></g>`;
    });
    s += `<text transform="translate(13 ${mt + ph / 2}) rotate(-90)" font-size="11" text-anchor="middle" fill="var(--faint)">Δh relative to mica (nm)</text>`;
    s += '</svg>';
    box.innerHTML = s;
  }

  // Box 2: schematic approach curves, rigid reference vs softer sample (not measured data)
  function stiffFigure() {
    const box = document.getElementById('useStiffFig');
    if (!box) return;
    const W = boxW(box), H = 250, ml = 50, mr = 14, mt = 14, mb = 40, pw = W - ml - mr, ph = H - mt - mb;
    const x0 = -20, x1 = 160, y0 = 0.965, y1 = 1.002;
    const X = (v) => ml + ((v - x0) / (x1 - x0)) * pw, Y = (v) => mt + (1 - (v - y0) / (y1 - y0)) * ph;
    // I/I_inf = 1 - 0.01 exp(-z/L), z from the 0.99 crossing: slope at the crossing = 0.01/L
    const curve = (L) => { const pts = []; for (let z = x0; z <= x1; z += 1) pts.push([z, 1 - 0.01 * Math.exp(-z / L)]); return pts; };
    const path = (pts) => pts.map(([x, y], i) => (i ? 'L' : 'M') + X(x).toFixed(1) + ' ' + Y(y).toFixed(1)).join('');
    const acc = 'var(--acc)';
    let s = `<svg viewBox="0 0 ${W} ${H}" xmlns="${NS}" font-family="Geist, Helvetica, Arial, sans-serif">`;
    s += `<defs><clipPath id="useStiffClip"><rect x="${ml}" y="${mt}" width="${pw}" height="${ph}"/></clipPath></defs>`;
    s += `<rect x="${ml}" y="${Y(0.99)}" width="${pw}" height="${Y(0.98) - Y(0.99)}" fill="var(--ov3)"/>`;
    s += `<text x="${ml + pw - 4}" y="${Y(0.98) - 4}" font-size="10.5" text-anchor="end" fill="var(--faint)">fit window 0.98–0.99</text>`;
    for (const t of [0.97, 0.98, 0.99, 1.0]) {
      s += `<line x1="${ml}" x2="${ml + pw}" y1="${Y(t)}" y2="${Y(t)}" stroke="var(--line-2)" stroke-width="1"/>`;
      s += `<text x="${ml - 6}" y="${Y(t) + 4}" font-size="11" text-anchor="end" fill="var(--faint)">${t.toFixed(2)}</text>`;
    }
    for (const t of [0, 50, 100, 150]) s += `<text x="${X(t)}" y="${mt + ph + 16}" font-size="11" text-anchor="middle" fill="var(--faint)">${t}</text>`;
    s += `<line x1="${ml}" x2="${ml}" y1="${mt}" y2="${mt + ph}" stroke="var(--muted)"/><line x1="${ml}" x2="${ml + pw}" y1="${mt + ph}" y2="${mt + ph}" stroke="var(--muted)"/>`;
    s += `<g clip-path="url(#useStiffClip)">`;
    s += `<path d="${path(curve(12))}" fill="none" stroke="var(--muted)" stroke-width="2" stroke-dasharray="6 4"><title>Rigid reference: steeper, slope s∞</title></path>`;
    s += `<path d="${path(curve(35))}" fill="none" stroke="${acc}" stroke-width="2"><title>Softer sample: less steep, slope s &lt; s∞</title></path>`;
    s += '</g>';
    s += `<text x="${X(-9)}" y="${Y(0.969)}" font-size="11.5" fill="var(--muted)">rigid reference, slope s∞ (dashed)</text>`;
    s += `<text x="${X(52)}" y="${Y(0.9955) + 14}" font-size="11.5" fill="${acc}" font-weight="600">softer sample, slope s</text>`;
    s += `<text x="${ml + pw / 2}" y="${H - 6}" font-size="11" text-anchor="middle" fill="var(--faint)">distance from the 0.99 crossing (nm)</text>`;
    s += `<text transform="translate(13 ${mt + ph / 2}) rotate(-90)" font-size="11" text-anchor="middle" fill="var(--faint)">I/I∞</text>`;
    s += '</svg>';
    box.innerHTML = s;
  }

  function init() {
    fillConventions(); chargeFigure(); stiffFigure();
    let t = null;
    window.addEventListener('resize', () => { clearTimeout(t); t = setTimeout(() => { chargeFigure(); stiffFigure(); }, 150); });
  }
  if (document.readyState === 'loading') document.addEventListener('DOMContentLoaded', init); else init();
})();
