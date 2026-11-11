/* CIHD - dual-bias SICM imaging module (SR-SICM Suite)
 * Leyla Nur Akay, Aarhus University, 2026.
 *
 * Processing follows the Python reference in examples/cihd/py (load.py, proc.py):
 *   1. height (nm) = raw int16 x data gain (Park header, tag 50435 @220) x 1000
 *   2. level():   plane fit + line-by-line median on the substrate, 4 iterations,
 *                 raised patches excluded by an iterative threshold, substrate = 0
 *   3. register:  cross-correlation of Gaussian-smoothed (sigma 2) patch masks,
 *                 10x upsampled DFT for the sub-pixel peak
 *                 (= skimage phase_cross_correlation, normalization=None, upsample_factor=10)
 *                 shift applied with bilinear interpolation, edges 'nearest'
 *                 (= scipy.ndimage.shift, order=1, mode='nearest')
 *   4. dh = h(-V) - h(+V), shown only where the two images overlap
 * Sign rule: dh > 0 -> more positive than the substrate (same rule as dz99).
 */
(function () {
  'use strict';

  const COL = { plus: '#d62728', minus: '#2f6fb5', dh: '#2ca02c' };
  // site root = folder above js/cihd.js (works from index.html and from tests/)
  const SITE_ROOT = (() => {
    const s = document.currentScript && document.currentScript.src;
    return s ? s.replace(/js\/cihd\.js(\?.*)?$/, '') : '';
  })();
  const EX_DIR = SITE_ROOT + 'examples/cihd/';

  /* =====================================================================
   * Numeric core (pure functions, exposed as CIHD.core for the tests)
   * ===================================================================== */

  function sortedCopy(a) { return Float64Array.from(a).sort(); }
  // numpy.percentile, default 'linear' method
  function pctSorted(s, p) {
    const n = s.length;
    if (!n) return NaN;
    const r = (p / 100) * (n - 1), i = Math.floor(r), f = r - i;
    return i + 1 < n ? s[i] + (s[i + 1] - s[i]) * f : s[i];
  }
  function percentile(a, p) { return pctSorted(sortedCopy(a), p); }
  function median(a) { return percentile(a, 50); }
  function stdPop(a) {           // numpy std, ddof = 0
    const n = a.length;
    if (!n) return NaN;
    let m = 0; for (let i = 0; i < n; i++) m += a[i]; m /= n;
    let s = 0; for (let i = 0; i < n; i++) s += (a[i] - m) ** 2;
    return Math.sqrt(s / n);
  }

  /* load.py level(): plane + line-by-line median leveling on the substrate;
   * iteratively excludes raised patches. Returns z (Float64Array) and the
   * substrate mask (Uint8Array, 1 = substrate). */
  function level(z, H, W, it = 4) {
    const N = H * W;
    let m = new Uint8Array(N).fill(1);
    const zp = new Float64Array(N);
    for (let k = 0; k < it; k++) {
      // plane fit (least squares on [x, y, 1]) over the current substrate mask
      let n = 0, mx = 0, my = 0, mz = 0;
      for (let y = 0, i = 0; y < H; y++) for (let x = 0; x < W; x++, i++)
        if (m[i]) { n++; mx += x; my += y; mz += z[i]; }
      mx /= n; my /= n; mz /= n;
      let sxx = 0, sxy = 0, syy = 0, sxz = 0, syz = 0;
      for (let y = 0, i = 0; y < H; y++) for (let x = 0; x < W; x++, i++)
        if (m[i]) {
          const dx = x - mx, dy = y - my, dz = z[i] - mz;
          sxx += dx * dx; sxy += dx * dy; syy += dy * dy; sxz += dx * dz; syz += dy * dz;
        }
      const det = sxx * syy - sxy * sxy;
      const a = det ? (sxz * syy - syz * sxy) / det : 0;
      const b = det ? (syz * sxx - sxz * sxy) / det : 0;
      const c = mz - a * mx - b * my;
      for (let y = 0, i = 0; y < H; y++) for (let x = 0; x < W; x++, i++) zp[i] = z[i] - (a * x + b * y + c);
      // line offsets (fast-scan rows) from substrate pixels only
      const row = new Float64Array(W);
      for (let y = 0; y < H; y++) {
        let c2 = 0;
        for (let x = 0; x < W; x++) if (m[y * W + x]) row[c2++] = zp[y * W + x];
        const off = c2 > 10 ? median(row.subarray(0, c2)) : 0;
        for (let x = 0; x < W; x++) zp[y * W + x] -= off;
      }
      // new substrate mask
      const s = sortedCopy(zp);
      const p20 = pctSorted(s, 20), p30 = pctSorted(s, 30), p60 = pctSorted(s, 60);
      const low = []; for (let i = 0; i < N; i++) if (zp[i] < p60) low.push(zp[i]);
      const thr = p30 + 3 * stdPop(low);
      const lim = Math.max(thr, p20 + 2.0);
      m = new Uint8Array(N);
      for (let i = 0; i < N; i++) m[i] = zp[i] < lim ? 1 : 0;
    }
    const sub = []; for (let i = 0; i < N; i++) if (m[i]) sub.push(zp[i]);
    const med = median(sub);
    for (let i = 0; i < N; i++) zp[i] -= med;
    return { z: zp, mask: m };
  }

  // scipy.ndimage.gaussian_filter, mode 'reflect', truncate 4.0
  function gaussian(img, H, W, sigma) {
    const r = Math.floor(4 * sigma + 0.5);
    const k = new Float64Array(2 * r + 1);
    let s = 0;
    for (let i = -r; i <= r; i++) { k[i + r] = Math.exp(-0.5 * i * i / (sigma * sigma)); s += k[i + r]; }
    for (let i = 0; i < k.length; i++) k[i] /= s;
    const refl = (i, n) => {
      if (n === 1) return 0;
      while (i < 0 || i >= n) { if (i < 0) i = -i - 1; if (i >= n) i = 2 * n - i - 1; }
      return i;
    };
    const tmp = new Float64Array(H * W), out = new Float64Array(H * W);
    for (let y = 0; y < H; y++) for (let x = 0; x < W; x++) {
      let v = 0; for (let j = -r; j <= r; j++) v += k[j + r] * img[y * W + refl(x + j, W)];
      tmp[y * W + x] = v;
    }
    for (let y = 0; y < H; y++) for (let x = 0; x < W; x++) {
      let v = 0; for (let j = -r; j <= r; j++) v += k[j + r] * tmp[refl(y + j, H) * W + x];
      out[y * W + x] = v;
    }
    return out;
  }

  // 1-D complex transform in place: radix-2 FFT for powers of two, direct DFT otherwise
  function transform1(re, im, inv) {
    const n = re.length, sg = inv ? 1 : -1;
    if (n & (n - 1)) {
      const or = new Float64Array(n), oi = new Float64Array(n);
      for (let k = 0; k < n; k++) {
        let sr = 0, si = 0;
        for (let t = 0; t < n; t++) {
          const a = sg * 2 * Math.PI * ((k * t) % n) / n, c = Math.cos(a), s = Math.sin(a);
          sr += re[t] * c - im[t] * s; si += re[t] * s + im[t] * c;
        }
        or[k] = sr; oi[k] = si;
      }
      re.set(or); im.set(oi);
      return;
    }
    for (let i = 1, j = 0; i < n; i++) {
      let bit = n >> 1;
      for (; j & bit; bit >>= 1) j ^= bit;
      j ^= bit;
      if (i < j) { let t = re[i]; re[i] = re[j]; re[j] = t; t = im[i]; im[i] = im[j]; im[j] = t; }
    }
    for (let len = 2; len <= n; len <<= 1) {
      const h = len >> 1;
      for (let kk = 0; kk < h; kk++) {
        const a = sg * 2 * Math.PI * kk / len, wr = Math.cos(a), wi = Math.sin(a);
        for (let i = 0; i < n; i += len) {
          const p = i + kk, q = p + h;
          const xr = re[q] * wr - im[q] * wi, xi = re[q] * wi + im[q] * wr;
          re[q] = re[p] - xr; im[q] = im[p] - xi; re[p] += xr; im[p] += xi;
        }
      }
    }
  }
  function transform2(re, im, H, W, inv) {
    const rr = new Float64Array(W), ri = new Float64Array(W);
    for (let y = 0; y < H; y++) {
      rr.set(re.subarray(y * W, y * W + W)); ri.set(im.subarray(y * W, y * W + W));
      transform1(rr, ri, inv);
      re.set(rr, y * W); im.set(ri, y * W);
    }
    const cr = new Float64Array(H), ci = new Float64Array(H);
    for (let x = 0; x < W; x++) {
      for (let y = 0; y < H; y++) { cr[y] = re[y * W + x]; ci[y] = im[y * W + x]; }
      transform1(cr, ci, inv);
      for (let y = 0; y < H; y++) { re[y * W + x] = cr[y]; im[y * W + x] = ci[y]; }
    }
    if (inv) { const s = 1 / (H * W); for (let i = 0; i < H * W; i++) { re[i] *= s; im[i] *= s; } }
  }
  const fftfreqInt = (n) => Array.from({ length: n }, (_, i) => (i < Math.ceil(n / 2) ? i : i - n));

  /* skimage.registration.phase_cross_correlation(ref, moving,
   *   upsample_factor=uf, normalization=None). Returns [dy, dx]: the shift to
   * apply to `moving` so that it lines up with `ref`. */
  function crossCorrShift(ref, mov, H, W, uf = 10) {
    const N = H * W;
    const ar = Float64Array.from(ref), ai = new Float64Array(N);
    const br = Float64Array.from(mov), bi = new Float64Array(N);
    transform2(ar, ai, H, W, false);
    transform2(br, bi, H, W, false);
    const pr = new Float64Array(N), pi = new Float64Array(N);      // ref * conj(mov)
    for (let i = 0; i < N; i++) { pr[i] = ar[i] * br[i] + ai[i] * bi[i]; pi[i] = ai[i] * br[i] - ar[i] * bi[i]; }
    const cr = Float64Array.from(pr), ci = Float64Array.from(pi);
    transform2(cr, ci, H, W, true);
    let best = -1, bi0 = 0;
    for (let i = 0; i < N; i++) { const v = cr[i] * cr[i] + ci[i] * ci[i]; if (v > best) { best = v; bi0 = i; } }
    let sy = Math.floor(bi0 / W), sx = bi0 % W;
    if (sy > Math.trunc(H / 2)) sy -= H;
    if (sx > Math.trunc(W / 2)) sx -= W;
    if (uf <= 1) return [sy, sx];
    // upsampled DFT around the coarse peak (skimage _upsampled_dft on conj(product))
    const ups = Math.ceil(uf * 1.5), dftshift = Math.trunc(ups / 2);
    const offY = dftshift - sy * uf, offX = dftshift - sx * uf;
    const fy = fftfreqInt(H), fx = fftfreqInt(W);
    // T[y][b] = sum_x conj(P[y,x]) * exp(-2 pi i (b - offX) fx[x] / (W uf))
    const Tr = new Float64Array(H * ups), Ti = new Float64Array(H * ups);
    const kc = new Float64Array(ups * W), ks = new Float64Array(ups * W);
    for (let b = 0; b < ups; b++) for (let x = 0; x < W; x++) {
      const a = -2 * Math.PI / (W * uf) * (b - offX) * fx[x];
      kc[b * W + x] = Math.cos(a); ks[b * W + x] = Math.sin(a);
    }
    for (let y = 0; y < H; y++) for (let b = 0; b < ups; b++) {
      let sr = 0, si = 0;
      for (let x = 0; x < W; x++) {
        const dr = pr[y * W + x], di = -pi[y * W + x], c = kc[b * W + x], s = ks[b * W + x];
        sr += dr * c - di * s; si += dr * s + di * c;
      }
      Tr[y * ups + b] = sr; Ti[y * ups + b] = si;
    }
    let bestU = -1, ba = 0, bb = 0;
    for (let a2 = 0; a2 < ups; a2++) {
      const rc = new Float64Array(H), rs = new Float64Array(H);
      for (let y = 0; y < H; y++) { const a = -2 * Math.PI / (H * uf) * (a2 - offY) * fy[y]; rc[y] = Math.cos(a); rs[y] = Math.sin(a); }
      for (let b = 0; b < ups; b++) {
        let sr = 0, si = 0;
        for (let y = 0; y < H; y++) {
          const tr = Tr[y * ups + b], ti = Ti[y * ups + b];
          sr += rc[y] * tr - rs[y] * ti; si += rc[y] * ti + rs[y] * tr;
        }
        const v = sr * sr + si * si;
        if (v > bestU) { bestU = v; ba = a2; bb = b; }
      }
    }
    return [sy + (ba - dftshift) / uf, sx + (bb - dftshift) / uf];
  }

  // scipy.ndimage.shift(img, (dy, dx), order=1, mode='nearest')
  function shiftBilinear(img, H, W, dy, dx) {
    const out = new Float64Array(H * W);
    for (let y = 0; y < H; y++) {
      let sy = Math.min(Math.max(y - dy, 0), H - 1);
      const y0 = Math.floor(sy), y1 = Math.min(y0 + 1, H - 1), fy = sy - y0;
      for (let x = 0; x < W; x++) {
        let sx = Math.min(Math.max(x - dx, 0), W - 1);
        const x0 = Math.floor(sx), x1 = Math.min(x0 + 1, W - 1), fx = sx - x0;
        out[y * W + x] = (1 - fy) * ((1 - fx) * img[y0 * W + x0] + fx * img[y0 * W + x1]) +
                         fy * ((1 - fx) * img[y1 * W + x0] + fx * img[y1 * W + x1]);
      }
    }
    return out;
  }

  // proc.py: valid region after the shift (rows y0..y1-1, cols x0..x1-1)
  function validRegion(H, W, dy, dx) {
    return { y0: Math.ceil(Math.max(dy, 0)), y1: H - Math.ceil(Math.max(-dy, 0)),
             x0: Math.ceil(Math.max(dx, 0)), x1: W - Math.ceil(Math.max(-dx, 0)) };
  }

  function crop(a, H, W, r) {
    const h = r.y1 - r.y0, w = r.x1 - r.x0, out = new a.constructor(h * w);
    for (let y = 0; y < h; y++) for (let x = 0; x < w; x++) out[y * w + x] = a[(y + r.y0) * W + x + r.x0];
    return out;
  }

  // numpy: rot90(fliplr(a) if flip else a, k)  (k quarter turns counter-clockwise)
  function orient(a, H, W, flip, k) {
    let cur = a, h = H, w = W;
    if (flip) { const o = new a.constructor(h * w); for (let y = 0; y < h; y++) for (let x = 0; x < w; x++) o[y * w + x] = cur[y * w + (w - 1 - x)]; cur = o; }
    for (let t = 0; t < ((k % 4) + 4) % 4; t++) {
      const o = new a.constructor(h * w), nh = w, nw = h;           // out[i][j] = in[j][w-1-i]
      for (let i = 0; i < nh; i++) for (let j = 0; j < nw; j++) o[i * nw + j] = cur[j * w + (w - 1 - i)];
      cur = o; h = nh; w = nw;
    }
    return { a: cur, H: h, W: w };
  }

  // scipy.ndimage.binary_erosion, cross structure, border_value 0
  function erode(mask, H, W, iter) {
    let cur = mask;
    for (let t = 0; t < iter; t++) {
      const o = new Uint8Array(H * W);
      for (let y = 1; y < H - 1; y++) for (let x = 1; x < W - 1; x++) {
        const i = y * W + x;
        o[i] = cur[i] && cur[i - 1] && cur[i + 1] && cur[i - W] && cur[i + W] ? 1 : 0;
      }
      cur = o;
    }
    return cur;
  }

  function maskedStats(D, mask) {
    const v = []; for (let i = 0; i < D.length; i++) if (mask[i]) v.push(D[i]);
    const s = sortedCopy(v);
    return { median: pctSorted(s, 50), q25: pctSorted(s, 25), q75: pctSorted(s, 75), n: s.length };
  }
  // figs.py rule: patch = P > p75(P), eroded twice
  function patchStatsP75(P, D, H, W) {
    const t = percentile(P, 75), m = new Uint8Array(H * W);
    for (let i = 0; i < m.length; i++) m[i] = P[i] > t ? 1 : 0;
    return maskedStats(D, erode(m, H, W, 2));
  }
  // levelling-mask rule: patch = not substrate (from level()), eroded twice
  function patchStatsMask(patch, D, H, W) { return maskedStats(D, erode(patch, H, W, 2)); }

  // mean of rows row-half .. row+half (clipped at the edges)
  function profileRows(a, H, W, row, half) {
    const y0 = Math.max(0, row - half), y1 = Math.min(H - 1, row + half), out = new Float64Array(W);
    for (let y = y0; y <= y1; y++) for (let x = 0; x < W; x++) out[x] += a[y * W + x];
    for (let x = 0; x < W; x++) out[x] /= (y1 - y0 + 1);
    return out;
  }

  const pyRound = (v) => {           // numpy/Python round half to even
    const f = Math.floor(v), d = v - f;
    return d > 0.5 ? f + 1 : d < 0.5 ? f : (f % 2 === 0 ? f : f + 1);
  };

  /* =====================================================================
   * File readers
   * ===================================================================== */

  /* Park XEI/XEP topography TIFF. Tag 50434 = raw data (int16, W x H),
   * tag 50435 = header: width int32 @100, height int32 @104, x scan size double @140 (um),
   * y @148, bias double @212 (V), data gain double @220. Height (nm) = raw x gain x 1000. */
  function readParkImage(buffer) {
    const dv = new DataView(buffer), n = buffer.byteLength;
    const bo = dv.getUint16(0);
    if (bo !== 0x4949 && bo !== 0x4D4D) throw new Error('not a TIFF file');
    const le = bo === 0x4949;
    if (dv.getUint16(2, le) !== 42) throw new Error('not a TIFF file');
    const ifd = dv.getUint32(4, le), cnt = dv.getUint16(ifd, le), tags = {};
    for (let i = 0; i < cnt; i++) {
      const o = ifd + 2 + 12 * i;
      if (o + 12 > n) break;
      tags[dv.getUint16(o, le)] = { type: dv.getUint16(o + 2, le), count: dv.getUint32(o + 4, le), value: dv.getUint32(o + 8, le) };
    }
    if (!tags[50434] || !tags[50435]) throw new Error('no Park image data (tags 50434/50435 missing). Approach-curve files belong in the curve modules.');
    const h = tags[50435].value;
    if (h + 228 > n) throw new Error('Park header truncated');
    const W = dv.getInt32(h + 100, le), H = dv.getInt32(h + 104, le);
    const xs = dv.getFloat64(h + 140, le), ys = dv.getFloat64(h + 148, le);
    const bias = dv.getFloat64(h + 212, le), gain = dv.getFloat64(h + 220, le);
    if (!(W > 0 && H > 0 && W * H <= 4096 * 4096)) throw new Error('image size in header not valid');
    const off = tags[50434].value, bytes = tags[50434].count;
    if (bytes < W * H * 2 || off + W * H * 2 > n) throw new Error('image data shorter than W x H int16');
    const z = new Float64Array(W * H), sc = gain * 1000;
    for (let i = 0; i < W * H; i++) z[i] = dv.getInt16(off + 2 * i, le) * sc;
    let date = '';
    const t = tags[306];
    if (t && t.count >= 19 && t.value + 19 <= n) for (let i = 0; i < 19; i++) date += String.fromCharCode(dv.getUint8(t.value + i));
    return { kind: 'park', z, W, H, xs, ys, bias, gain, date };
  }

  /* plain text / CSV matrix: one image row per line; comma, semicolon, tab or space separated.
   * Lines that are not purely numeric (headers, comments) are skipped. Values in nm. */
  function readTextMatrix(text) {
    const rows = [];
    for (const line of text.split(/\r?\n/)) {
      const s = line.trim();
      if (!s) continue;
      const tok = s.split(/[,;\t ]+/);
      const v = tok.map(Number);
      if (v.some((x) => !isFinite(x))) continue;
      rows.push(v);
    }
    if (rows.length < 4) throw new Error('fewer than 4 numeric rows');
    const W = rows[0].length;
    const good = rows.filter((r) => r.length === W);
    if (good.length !== rows.length) throw new Error('rows have different lengths');
    const H = good.length, z = new Float64Array(W * H);
    good.forEach((r, y) => r.forEach((v, x) => { z[y * W + x] = v; }));
    return { kind: 'text', z, W, H, xs: NaN, ys: NaN, bias: NaN, gain: NaN, date: '' };
  }

  /* =====================================================================
   * Drawing helpers
   * ===================================================================== */

  const clamp01 = (v) => (v < 0 ? 0 : v > 1 ? 1 : v);
  const CMAPS = {
    afmhot: (t) => [clamp01(2 * t), clamp01(2 * t - 0.5), clamp01(2 * t - 1)],
    gray: (t) => [t, t, t],
  };
  function mapImage(a, H, W, vmin, vmax, cmap) {
    const cv = document.createElement('canvas'); cv.width = W; cv.height = H;
    const ctx = cv.getContext('2d'), im = ctx.createImageData(W, H), f = CMAPS[cmap];
    const span = vmax - vmin || 1;
    for (let i = 0; i < W * H; i++) {
      const c = f(clamp01((a[i] - vmin) / span));
      im.data[4 * i] = 255 * c[0]; im.data[4 * i + 1] = 255 * c[1]; im.data[4 * i + 2] = 255 * c[2]; im.data[4 * i + 3] = 255;
    }
    ctx.putImageData(im, 0, 0);
    return cv.toDataURL('image/png');
  }
  function cbarImage(cmap) {
    const cv = document.createElement('canvas'); cv.width = 1; cv.height = 256;
    const ctx = cv.getContext('2d'), im = ctx.createImageData(1, 256), f = CMAPS[cmap];
    for (let i = 0; i < 256; i++) { const c = f(1 - i / 255); im.data[4 * i] = 255 * c[0]; im.data[4 * i + 1] = 255 * c[1]; im.data[4 * i + 2] = 255 * c[2]; im.data[4 * i + 3] = 255; }
    ctx.putImageData(im, 0, 0);
    return cv.toDataURL('image/png');
  }
  function niceTicks(lo, hi, want) {
    if (!(hi > lo)) return [lo];
    const raw = (hi - lo) / Math.max(1, want), mag = Math.pow(10, Math.floor(Math.log10(raw)));
    const step = [1, 2, 2.5, 5, 10].map((m) => m * mag).find((s) => s >= raw) || 10 * mag;
    const out = [];
    for (let v = Math.ceil(lo / step - 1e-9) * step; v <= hi + 1e-9 * step; v += step) out.push(Math.abs(v) < step * 1e-9 ? 0 : v);
    return out;
  }
  const fmt = (v, d) => (isFinite(v) ? v.toFixed(d) : '—');
  const fmtTick = (v) => { const r = Math.round(v * 1000) / 1000; return String(r).replace('-', '−'); };
  const esc = (s) => String(s).replace(/[&<>"']/g, (c) => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' }[c]));
  const signed = (v, d) => (isFinite(v) ? (v > 0 ? '+' : v < 0 ? '−' : '') + Math.abs(v).toFixed(d) : '—');

  function themeColors(exportMode) {
    if (exportMode) return { bg: '#ffffff', txt: '#1a1a1a', muted: '#555555', line: '#bbbbbb', grid: '#e6e6e6' };
    const cs = getComputedStyle(document.documentElement);
    const g = (n, d) => (cs.getPropertyValue(n) || '').trim() || d;
    return { bg: 'none', txt: g('--txt', '#e8e8e8'), muted: g('--muted', '#9aa0ab'), line: g('--hair', 'rgba(255,255,255,.2)'), grid: g('--line-2', 'rgba(255,255,255,.06)') };
  }

  /* one map panel as an SVG string. Geometry in user units; the SVG scales with its container. */
  function mapGeom(H, W, px, py) {
    const ML = 40, MT = 24, MB = 34, CBG = 8, CBW = 9, MR = 44;
    let PW = 270, PH = PW * (H * py) / (W * px);
    if (PH > 300) { PH = 300; PW = PH * (W * px) / (H * py); }
    return { ML, MT, MB, CBG, CBW, MR, PW, PH, VW: ML + PW + CBG + CBW + MR, VH: MT + PH + MB };
  }
  function mapSVG(o, th) {
    const { H, W, px, py } = o, g = mapGeom(H, W, px, py);
    const wum = W * px, hum = H * py;
    const X = (um) => g.ML + (um / wum) * g.PW, Y = (um) => g.MT + (um / hum) * g.PH;
    const fs = 'font-family="Geist, Helvetica, Arial, sans-serif"';
    let s = `<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${g.VW.toFixed(1)} ${g.VH.toFixed(1)}" ${fs} class="cihd-map-svg" role="img" aria-label="${esc(o.title)} map">`;
    if (th.bg !== 'none') s += `<rect width="100%" height="100%" fill="${th.bg}"/>`;
    s += `<text x="${g.ML}" y="15" font-size="11.5" font-weight="600" fill="${th.txt}">${esc(o.title)}</text>`;
    s += `<image href="${o.img}" x="${g.ML}" y="${g.MT}" width="${g.PW.toFixed(2)}" height="${g.PH.toFixed(2)}" preserveAspectRatio="none" style="image-rendering:pixelated"/>`;
    s += `<rect x="${g.ML}" y="${g.MT}" width="${g.PW.toFixed(2)}" height="${g.PH.toFixed(2)}" fill="none" stroke="${th.line}" stroke-width="0.8"/>`;
    // badge
    if (o.badge) {
      const bw = 8 + 6.2 * o.badge.length;
      s += `<rect x="${(g.ML + g.PW - bw - 4).toFixed(1)}" y="${g.MT + 4}" width="${bw.toFixed(1)}" height="15" rx="2" fill="${o.badgeBg}" fill-opacity="0.9"/>`;
      s += `<text x="${(g.ML + g.PW - 8).toFixed(1)}" y="${g.MT + 15}" font-size="10" font-weight="600" text-anchor="end" fill="${o.badgeFg}">${esc(o.badge)}</text>`;
    }
    // profile band + line
    const rowY = (r) => g.MT + ((r + 0.5) / H) * g.PH;
    const y0 = Math.max(0, o.row - o.half), y1 = Math.min(H - 1, o.row + o.half);
    s += `<rect class="cihd-band" x="${g.ML}" y="${(g.MT + (y0 / H) * g.PH).toFixed(2)}" width="${g.PW.toFixed(2)}" height="${(((y1 - y0 + 1) / H) * g.PH).toFixed(2)}" fill="${o.lineColor}" fill-opacity="0.18"/>`;
    s += `<line class="cihd-line" x1="${g.ML}" x2="${(g.ML + g.PW).toFixed(2)}" y1="${rowY(o.row).toFixed(2)}" y2="${rowY(o.row).toFixed(2)}" stroke="${o.lineColor}" stroke-width="1.6" stroke-dasharray="5 3"/>`;
    if (o.cropRect) {
      const c = o.cropRect;
      s += `<rect class="cihd-croprect" x="${(g.ML + (c.x0 / W) * g.PW).toFixed(2)}" y="${(g.MT + (c.y0 / H) * g.PH).toFixed(2)}" width="${(((c.x1 - c.x0) / W) * g.PW).toFixed(2)}" height="${(((c.y1 - c.y0) / H) * g.PH).toFixed(2)}" fill="none" stroke="#ffffff" stroke-width="1.2" stroke-dasharray="3 2"/>`;
    }
    // axes ticks (um)
    for (const t of niceTicks(0, wum, 4)) {
      const x = X(t).toFixed(2);
      s += `<line x1="${x}" x2="${x}" y1="${(g.MT + g.PH).toFixed(2)}" y2="${(g.MT + g.PH + 3).toFixed(2)}" stroke="${th.muted}" stroke-width="0.8"/>`;
      s += `<text x="${x}" y="${(g.MT + g.PH + 13).toFixed(2)}" font-size="9.5" text-anchor="middle" fill="${th.muted}">${fmtTick(t)}</text>`;
    }
    for (const t of niceTicks(0, hum, 3)) {
      const y = Y(t).toFixed(2);
      s += `<line x1="${g.ML - 3}" x2="${g.ML}" y1="${y}" y2="${y}" stroke="${th.muted}" stroke-width="0.8"/>`;
      s += `<text x="${g.ML - 5}" y="${(Y(t) + 3.2).toFixed(2)}" font-size="9.5" text-anchor="end" fill="${th.muted}">${fmtTick(t)}</text>`;
    }
    s += `<text x="${(g.ML + g.PW / 2).toFixed(2)}" y="${(g.VH - 4).toFixed(2)}" font-size="10" text-anchor="middle" fill="${th.muted}">x (µm)</text>`;
    s += `<text transform="translate(11 ${(g.MT + g.PH / 2).toFixed(2)}) rotate(-90)" font-size="10" text-anchor="middle" fill="${th.muted}">y (µm)</text>`;
    // colour bar
    const cx = g.ML + g.PW + g.CBG;
    s += `<image href="${o.cbar}" x="${cx.toFixed(2)}" y="${g.MT}" width="${g.CBW}" height="${g.PH.toFixed(2)}" preserveAspectRatio="none"/>`;
    s += `<rect x="${cx.toFixed(2)}" y="${g.MT}" width="${g.CBW}" height="${g.PH.toFixed(2)}" fill="none" stroke="${th.line}" stroke-width="0.6"/>`;
    for (const t of niceTicks(o.vmin, o.vmax, 4)) {
      const y = g.MT + (1 - (t - o.vmin) / (o.vmax - o.vmin)) * g.PH;
      s += `<line x1="${(cx + g.CBW).toFixed(2)}" x2="${(cx + g.CBW + 3).toFixed(2)}" y1="${y.toFixed(2)}" y2="${y.toFixed(2)}" stroke="${th.muted}" stroke-width="0.8"/>`;
      s += `<text x="${(cx + g.CBW + 5).toFixed(2)}" y="${(y + 3.2).toFixed(2)}" font-size="9.5" fill="${th.muted}">${fmtTick(t)}</text>`;
    }
    s += `<text transform="translate(${(g.VW - 4).toFixed(2)} ${(g.MT + g.PH / 2).toFixed(2)}) rotate(-90)" font-size="10" text-anchor="middle" fill="${th.muted}">${esc(o.unit)}</text>`;
    s += '</svg>';
    return s;
  }

  function profileSVG(p, width, th, opt = {}) {
    const VW = Math.max(300, Math.round(width)), VH = opt.height || 270;
    const ML = 50, MR = 14, MT = 30, MB = 40, PW = VW - ML - MR, PH = VH - MT - MB;
    const all = [...p.plus, ...p.minus, ...p.dh, 0];
    let lo = Math.min(...all), hi = Math.max(...all);
    const pad = 0.06 * (hi - lo || 1); lo -= pad; hi += pad;
    const x1 = p.x[p.x.length - 1] || 1;
    const X = (v) => ML + (v / x1) * PW, Y = (v) => MT + (1 - (v - lo) / (hi - lo)) * PH;
    const fs = 'font-family="Geist, Helvetica, Arial, sans-serif"';
    let s = `<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${VW} ${VH}" ${fs} class="cihd-prof-svg" role="img" aria-label="Line profile">`;
    if (th.bg !== 'none') s += `<rect width="100%" height="100%" fill="${th.bg}"/>`;
    for (const t of niceTicks(lo, hi, 5)) {
      const y = Y(t).toFixed(2);
      s += `<line x1="${ML}" x2="${ML + PW}" y1="${y}" y2="${y}" stroke="${th.grid}" stroke-width="0.8"/>`;
      s += `<text x="${ML - 6}" y="${(Y(t) + 3.5).toFixed(2)}" font-size="10.5" text-anchor="end" fill="${th.muted}">${fmtTick(t)}</text>`;
    }
    for (const t of niceTicks(0, x1, 8)) {
      const x = X(t).toFixed(2);
      s += `<line x1="${x}" x2="${x}" y1="${MT + PH}" y2="${MT + PH + 4}" stroke="${th.muted}" stroke-width="0.8"/>`;
      s += `<text x="${x}" y="${MT + PH + 16}" font-size="10.5" text-anchor="middle" fill="${th.muted}">${fmtTick(t)}</text>`;
    }
    s += `<line x1="${ML}" x2="${ML + PW}" y1="${Y(0).toFixed(2)}" y2="${Y(0).toFixed(2)}" stroke="#999999" stroke-width="0.9"/>`;
    s += `<line x1="${ML}" x2="${ML}" y1="${MT}" y2="${MT + PH}" stroke="${th.muted}" stroke-width="0.8"/>`;
    s += `<line x1="${ML}" x2="${ML + PW}" y1="${MT + PH}" y2="${MT + PH}" stroke="${th.muted}" stroke-width="0.8"/>`;
    const path = (arr) => arr.map((v, i) => (i ? 'L' : 'M') + X(p.x[i]).toFixed(2) + ' ' + Y(v).toFixed(2)).join('');
    s += `<path d="${path(p.plus)}" fill="none" stroke="${COL.plus}" stroke-width="1.6" stroke-linejoin="round"/>`;
    s += `<path d="${path(p.minus)}" fill="none" stroke="${COL.minus}" stroke-width="1.6" stroke-linejoin="round"/>`;
    // dh dashed: red and green are hard to tell apart with deuteranopia, so the line style differs too
    s += `<path d="${path(p.dh)}" fill="none" stroke="${COL.dh}" stroke-width="1.8" stroke-dasharray="6 3" stroke-linejoin="round"/>`;
    // legend
    let lx = ML + 4;
    for (const [lab, c, dash] of [[p.labels.plus, COL.plus, ''], [p.labels.minus, COL.minus, ''], ['Δh', COL.dh, '6 3']]) {
      s += `<line x1="${lx}" x2="${lx + 18}" y1="14" y2="14" stroke="${c}" stroke-width="2"${dash ? ` stroke-dasharray="${dash}"` : ''}/>`;
      s += `<text x="${lx + 23}" y="18" font-size="11" fill="${th.txt}">${esc(lab)}</text>`;
      lx += 34 + 6.3 * lab.length;
    }
    s += `<text x="${ML + PW / 2}" y="${VH - 6}" font-size="11" text-anchor="middle" fill="${th.muted}">Position along the line (µm)</text>`;
    s += `<text transform="translate(13 ${MT + PH / 2}) rotate(-90)" font-size="11" text-anchor="middle" fill="${th.muted}">Height (nm)</text>`;
    s += `<g class="cihd-hover" style="display:none"><line y1="${MT}" y2="${MT + PH}" stroke="${th.muted}" stroke-width="0.8" stroke-dasharray="2 2"/><text font-size="10.5" fill="${th.txt}"></text></g>`;
    s += `<rect class="cihd-hit" x="${ML}" y="${MT}" width="${PW}" height="${PH}" fill="transparent"/>`;
    s += '</svg>';
    return { svg: s, geom: { ML, PW, x1 } };
  }

  /* =====================================================================
   * Module state and UI
   * ===================================================================== */

  const S = {
    mode: 'example',
    manifest: null,
    exId: 'DPTAP',
    exCache: {},
    user: { files: [], pipNeg: false, iters: 4, base: null, manual: null, flip: false, rot: 0, crop: null, cropMode: false },
    view: null,
  };
  let root = null;
  const $ = (sel) => root.querySelector(sel);

  function msg(text, kind) {
    const el = $('#cihdMsg');
    if (!el) return;
    el.textContent = text || '';
    el.className = 'cihd-msg' + (kind ? ' cihd-msg-' + kind : '');
  }

  function buildUI() {
    root.innerHTML = `
      <div class="cihd-head">
        <div>
          <h3 class="cihd-title">Dual-bias imaging (CIHD)</h3>
          <p class="cihd-sub">Topography recorded at +V and −V over the same area. The height difference
            Δh = h(−V) − h(+V) has the sign of the surface charge relative to the substrate (same rule as Δz99).</p>
        </div>
        <div class="cihd-seg" role="group" aria-label="Data source">
          <button type="button" data-mode="example" aria-pressed="true">Example</button>
          <button type="button" data-mode="user" aria-pressed="false">Your data</button>
        </div>
      </div>

      <div class="cihd-permission" id="cihdPermission" hidden></div>

      <div class="cihd-panel" id="cihdExample">
        <div class="cihd-row">
          <div class="cihd-seg" id="cihdLipids" role="group" aria-label="Example"></div>
          <label class="cihd-check" title="In these 2015 Park XEP files the paper's +100 mV is header bias −0.10 V">
            <input type="checkbox" checked disabled> Pipette = −header bias (2015 data)
          </label>
        </div>
        <p class="cihd-note" id="cihdExNote"></p>
      </div>

      <div class="cihd-panel" id="cihdUser" hidden>
        <div class="cihd-drop" id="cihdDrop" tabindex="0" role="button" aria-label="Load two topography images">
          <input type="file" id="cihdFile" multiple accept=".tif,.tiff,.txt,.csv,.dat,.asc" hidden>
          <strong>Load two topography images of the same area</strong> (one at +V, one at −V)<br>
          <span>Park XEI/XEP .tiff, or a text/CSV height matrix in nm. Drop files here or click.</span>
        </div>
        <div class="cihd-files" id="cihdFiles"></div>
        <div class="cihd-row">
          <label class="cihd-check"><input type="checkbox" id="cihdPipNeg"> Pipette = −header bias (2015 data)</label>
          <span class="cihd-hint" id="cihdPipHint"></span>
        </div>
        <div id="cihdSteps" hidden>
          <details class="cihd-step" open>
            <summary>a) Levelling</summary>
            <p>Plane fit and line-by-line median on the substrate; raised patches are excluded by an iterative
              threshold; substrate = 0.</p>
            <div class="cihd-row">
              <label class="cihd-inline">Iterations <input type="number" id="cihdIters" min="1" max="10" step="1" value="4"></label>
              <span class="cihd-readout" id="cihdLevelOut"></span>
            </div>
          </details>
          <details class="cihd-step" open>
            <summary>b) Registration</summary>
            <p>The −V image is shifted onto the +V image. Shift from the cross-correlation of the patch masks
              (Gaussian σ = 2 px), 0.1 px resolution. Fine-tune by hand if needed.</p>
            <div class="cihd-row">
              <span class="cihd-readout" id="cihdShiftAuto"></span>
            </div>
            <div class="cihd-row">
              <label class="cihd-inline">Δy (px) <input type="number" id="cihdDy" step="0.1"></label>
              <label class="cihd-inline">Δx (px) <input type="number" id="cihdDx" step="0.1"></label>
              <button type="button" class="btn" id="cihdShiftReset">Use automatic</button>
            </div>
          </details>
          <details class="cihd-step">
            <summary>c) Δh and overlap</summary>
            <p>Δh = h(−V) − h(+V), computed only where the two images overlap after the shift.</p>
            <span class="cihd-readout" id="cihdOverlap"></span>
          </details>
          <details class="cihd-step">
            <summary>d) Orientation and crop (optional)</summary>
            <div class="cihd-row">
              <label class="cihd-check"><input type="checkbox" id="cihdFlip"> Mirror left–right</label>
              <button type="button" class="btn" id="cihdRot">Rotate 90°</button>
              <span class="cihd-readout" id="cihdRotOut">0°</span>
            </div>
            <div class="cihd-row">
              <button type="button" class="btn" id="cihdCropMode" aria-pressed="false">Draw crop box on a map</button>
              <button type="button" class="btn" id="cihdCropReset">Reset crop</button>
              <span class="cihd-readout" id="cihdCropOut"></span>
            </div>
          </details>
        </div>
      </div>

      <div class="cihd-msg" id="cihdMsg" role="status" aria-live="polite"></div>

      <div id="cihdView" hidden>
        <div class="cihd-maps">
          <figure class="cihd-map" data-key="P"></figure>
          <figure class="cihd-map" data-key="M"></figure>
          <figure class="cihd-map" data-key="D"></figure>
        </div>
        <div class="cihd-row cihd-profctl">
          <label class="cihd-inline cihd-grow">Line position
            <input type="range" id="cihdRow" min="0" max="1" step="1" value="0">
            <span class="cihd-readout" id="cihdRowOut"></span></label>
          <label class="cihd-inline">Width
            <select id="cihdHalf">
              <option value="0">1 line</option><option value="1">3 lines</option>
              <option value="2" selected>5 lines (±2)</option><option value="3">7 lines</option><option value="4">9 lines</option>
            </select></label>
        </div>
        <div class="cihd-prof" id="cihdProf"></div>
        <p class="cihd-caption" id="cihdCaption"></p>
        <div class="cihd-summary" id="cihdSummary"></div>
        <div class="cihd-row cihd-export">
          <span class="cihd-label">Export</span>
          <button type="button" class="btn" data-exp="maps">Maps (CSV)</button>
          <button type="button" class="btn" data-exp="profile">Profile (CSV)</button>
          <button type="button" class="btn" data-exp="png">Figure (PNG)</button>
          <button type="button" class="btn" data-exp="svg">Figure (SVG)</button>
        </div>
      </div>`;

    root.querySelectorAll('[data-mode]').forEach((b) => b.addEventListener('click', () => setMode(b.dataset.mode)));
    const drop = $('#cihdDrop'), inp = $('#cihdFile');
    drop.addEventListener('click', () => inp.click());
    drop.addEventListener('keydown', (e) => { if (e.key === 'Enter' || e.key === ' ') { e.preventDefault(); inp.click(); } });
    drop.addEventListener('dragover', (e) => { e.preventDefault(); drop.classList.add('is-over'); });
    drop.addEventListener('dragleave', () => drop.classList.remove('is-over'));
    drop.addEventListener('drop', (e) => { e.preventDefault(); drop.classList.remove('is-over'); loadUserFiles(e.dataTransfer.files); });
    inp.addEventListener('change', () => { loadUserFiles(inp.files); inp.value = ''; });
    $('#cihdPipNeg').addEventListener('change', (e) => { S.user.pipNeg = e.target.checked; autoRoles(); renderFiles(); runUser(); });
    $('#cihdIters').addEventListener('change', (e) => { S.user.iters = Math.max(1, Math.min(10, parseInt(e.target.value, 10) || 4)); runUser(); });
    const onShift = () => {
      const dy = parseFloat($('#cihdDy').value), dx = parseFloat($('#cihdDx').value);
      if (isFinite(dy) && isFinite(dx)) { S.user.manual = [dy, dx]; S.user.crop = null; composeUser(); }
    };
    $('#cihdDy').addEventListener('change', onShift);
    $('#cihdDx').addEventListener('change', onShift);
    $('#cihdShiftReset').addEventListener('click', () => { S.user.manual = null; S.user.crop = null; composeUser(); });
    $('#cihdFlip').addEventListener('change', (e) => { S.user.flip = e.target.checked; S.user.crop = null; composeUser(); });
    $('#cihdRot').addEventListener('click', () => { S.user.rot = (S.user.rot + 1) % 4; S.user.crop = null; composeUser(); });
    $('#cihdCropMode').addEventListener('click', () => {
      S.user.cropMode = !S.user.cropMode;
      $('#cihdCropMode').setAttribute('aria-pressed', String(S.user.cropMode));
      root.classList.toggle('cihd-cropping', S.user.cropMode);
    });
    $('#cihdCropReset').addEventListener('click', () => { S.user.crop = null; composeUser(); });
    $('#cihdRow').addEventListener('input', (e) => setRow(parseInt(e.target.value, 10)));
    $('#cihdHalf').addEventListener('change', (e) => { if (S.view) { S.view.half = parseInt(e.target.value, 10); updateLine(); } });
    root.querySelectorAll('[data-exp]').forEach((b) => b.addEventListener('click', () => doExport(b.dataset.exp)));
    if (window.ResizeObserver) new ResizeObserver(() => { if (S.view && !root.hidden) renderProfile(); }).observe($('#cihdProf'));
    window.addEventListener('sicm-theme', () => { if (S.view) { renderMaps(); renderProfile(); } });
  }

  function setMode(mode) {
    S.mode = mode;
    root.querySelectorAll('[data-mode]').forEach((b) => b.setAttribute('aria-pressed', String(b.dataset.mode === mode)));
    $('#cihdExample').hidden = mode !== 'example';
    $('#cihdUser').hidden = mode !== 'user';
    showPermission();
    msg('');
    if (mode === 'example') loadExample(S.exId);
    else if (S.user.base) composeUser();
    else if (S.user.files.length) runUser();
    else { S.view = null; $('#cihdView').hidden = true; }
  }

  function showPermission() {
    const el = $('#cihdPermission'), p = S.manifest && S.manifest.permission;
    const show = S.mode === 'example' && p && p.status !== 'public';
    el.hidden = !show;
    if (show) el.textContent = p.note || 'Preliminary data.';
  }

  /* ---------------- example mode ---------------- */
  async function fetchText(url) {
    const r = await fetch(url, { cache: 'no-cache' });
    if (!r.ok) throw new Error(url + ': HTTP ' + r.status);
    return r.text();
  }
  async function ensureManifest() {
    if (S.manifest) return S.manifest;
    S.manifest = JSON.parse(await fetchText(EX_DIR + 'manifest.json'));
    const box = $('#cihdLipids');
    box.innerHTML = S.manifest.lipids.map((l) => `<button type="button" data-ex="${esc(l.id)}" aria-pressed="false">${esc(l.id)} <span>${esc(l.charge)}</span></button>`).join('');
    box.querySelectorAll('[data-ex]').forEach((b) => b.addEventListener('click', () => loadExample(b.dataset.ex)));
    return S.manifest;
  }

  async function loadExample(id) {
    try {
      msg('Loading example …');
      const man = await ensureManifest();
      showPermission();
      const L = man.lipids.find((l) => l.id === id) || man.lipids[0];
      S.exId = L.id;
      root.querySelectorAll('[data-ex]').forEach((b) => b.setAttribute('aria-pressed', String(b.dataset.ex === L.id)));
      if (!S.exCache[L.id]) {
        const get = async (k) => readTextMatrix(await fetchText(`${EX_DIR}CIHD_${L.id}_map_${k}_nm.csv`));
        const [P, M, D] = await Promise.all([get('plus100mV'), get('minus100mV'), get('dh')]);
        if (P.W !== M.W || P.W !== D.W || P.H !== M.H || P.H !== D.H) throw new Error('example maps differ in size');
        const lv = level(P.z, P.H, P.W, 4);
        const patch = lv.mask.map((v) => 1 - v);
        S.exCache[L.id] = { P: P.z, M: M.z, D: D.z, H: P.H, W: P.W, patch };
      }
      const c = S.exCache[L.id], px = L.px_nm / 1000;
      const pub = L.published_dh_nm;
      $('#cihdExNote').innerHTML =
        `<strong>${esc(L.id)}</strong> (${esc(L.charge)}) patches on mica, imaged at +100 mV and −100 mV ` +
        `(${(c.W * px).toFixed(1)} × ${(c.H * px).toFixed(1)} µm, ${L.px_nm.toFixed(1)} nm/px). ` +
        `Expected: Δh ${pub.median > 0 ? '> 0, patches more positive than mica' : '< 0, patches more negative than mica'}; ` +
        `published median ${signed(pub.median, 1)} nm.`;
      setView({
        name: L.id, source: 'example', P: c.P, M: c.M, D: c.D, H: c.H, W: c.W, px, py: px, patch: c.patch,
        labels: { plus: '+100 mV', minus: '−100 mV', plusCsv: 'plus100mV', minusCsv: 'minus100mV' },
        row: Math.min(c.H - 1, pyRound(L.row_um / px)), half: man.profile_rows ? (man.profile_rows - 1) / 2 : 2,
        published: pub,
        caption: `${man.source} Δh = h(−100 mV) − h(+100 mV), relative to the mica substrate. ` +
                 `Profile: mean of the selected lines. Bias labels follow the paper (pipette potential = −header bias in these 2015 files).`,
      });
      msg('');
    } catch (e) {
      const local = location.protocol === 'file:';
      msg((local ? 'Examples need a web server (open the site via http://, not file://). ' : 'Could not load the example: ') + e.message, 'err');
    }
  }

  /* ---------------- user data mode ---------------- */
  async function loadUserFiles(fileList) {
    const files = Array.from(fileList || []);
    if (!files.length) return;
    const read = [];
    for (const f of files) {
      try {
        const isTif = /\.tiff?$/i.test(f.name);
        const d = isTif ? readParkImage(await f.arrayBuffer()) : readTextMatrix(await f.text());
        d.name = f.name;
        d.pxNm = isTif ? (d.xs / d.W) * 1000 : NaN;
        d.pyNm = isTif ? (d.ys / d.H) * 1000 : NaN;
        read.push(d);
      } catch (e) {
        msg(`${f.name}: ${e.message}`, 'err');
      }
    }
    if (!read.length) return;
    S.user.files = S.user.files.concat(read).slice(-2);
    S.user.base = null; S.user.manual = null; S.user.crop = null;
    const old = S.user.files.some((f) => /^(19|20)(0\d|1[0-5])/.test(f.date || ''));
    $('#cihdPipHint').textContent = old ? 'Files from 2015 or earlier: in the Klausen et al. 2016 data the paper\'s +V is header −V; tick this option for those files.' : '';
    autoRoles();
    renderFiles();
    runUser();
  }

  // role from the header bias (sign x convention); text files: first = +V, second = −V
  function autoRoles() {
    const F = S.user.files;
    F.forEach((f, i) => {
      if (isFinite(f.bias) && f.bias !== 0) {
        const v = S.user.pipNeg ? -f.bias : f.bias;
        f.role = v > 0 ? 'plus' : 'minus';
      } else if (!f.role) f.role = i === 0 ? 'plus' : 'minus';
    });
  }

  function renderFiles() {
    const box = $('#cihdFiles');
    box.innerHTML = S.user.files.map((f, i) => {
      const hdr = isFinite(f.bias) ? `header bias ${f.bias >= 0 ? '+' : '−'}${Math.abs(f.bias).toFixed(3)} V` : 'no bias in file';
      const size = f.kind === 'park' ? `${f.W} × ${f.H} px · ${f.xs.toFixed(2)} × ${f.ys.toFixed(2)} µm · ${f.pxNm.toFixed(1)} nm/px` : `${f.W} × ${f.H} px`;
      const pxIn = f.kind === 'text'
        ? `<label class="cihd-inline">Pixel size <input type="number" step="0.1" min="0" data-px="${i}" value="${isFinite(f.pxNm) ? f.pxNm : ''}" placeholder="nm"> nm</label>` : '';
      return `<div class="cihd-file">
          <div class="cihd-file-name" title="${esc(f.name)}">${esc(f.name)}</div>
          <div class="cihd-file-meta">${size} · ${hdr}${f.date ? ' · ' + esc(f.date) : ''}</div>
          <div class="cihd-row">
            <div class="cihd-seg cihd-seg-sm" role="group" aria-label="Role of ${esc(f.name)}">
              <button type="button" data-role="plus" data-i="${i}" aria-pressed="${f.role === 'plus'}" class="is-plus">+V</button>
              <button type="button" data-role="minus" data-i="${i}" aria-pressed="${f.role === 'minus'}" class="is-minus">−V</button>
            </div>
            ${pxIn}
            <button type="button" class="btn cihd-x" data-rm="${i}" aria-label="Remove ${esc(f.name)}">Remove</button>
          </div>
        </div>`;
    }).join('');
    box.querySelectorAll('[data-role]').forEach((b) => b.addEventListener('click', () => {
      const i = +b.dataset.i, F = S.user.files;
      F[i].role = b.dataset.role;
      if (F.length === 2) F[1 - i].role = b.dataset.role === 'plus' ? 'minus' : 'plus';
      S.user.base = null; renderFiles(); runUser();
    }));
    box.querySelectorAll('[data-px]').forEach((inp) => inp.addEventListener('change', () => {
      const v = parseFloat(inp.value), F = S.user.files;
      F.forEach((f) => { if (f.kind === 'text') { f.pxNm = v; f.pyNm = v; } });   // same pixel size for both
      S.user.base = null; renderFiles(); runUser();
    }));
    box.querySelectorAll('[data-rm]').forEach((b) => b.addEventListener('click', () => {
      S.user.files.splice(+b.dataset.rm, 1); S.user.base = null;
      renderFiles(); runUser();
    }));
  }

  function runUser() {
    const F = S.user.files;
    $('#cihdSteps').hidden = true;
    if (S.mode !== 'user') return;
    if (F.length < 2) { $('#cihdView').hidden = true; S.view = null; msg(F.length ? 'Load the second image.' : ''); return; }
    const fp = F.find((f) => f.role === 'plus'), fm = F.find((f) => f.role === 'minus');
    if (!fp || !fm) { msg('Assign one image to +V and the other to −V.', 'err'); return; }
    if (fp.W !== fm.W || fp.H !== fm.H) { msg(`The two images differ in size (${fp.W}×${fp.H} vs ${fm.W}×${fm.H} px); they must cover the same area with the same pixels.`, 'err'); return; }
    if (!(fp.pxNm > 0)) { msg('Enter the pixel size (nm) for the text files.', 'err'); return; }
    if (fp.kind === 'park' && fm.kind === 'park' && (Math.abs(fp.xs - fm.xs) > 1e-6 || Math.abs(fp.ys - fm.ys) > 1e-6)) {
      msg('The two images have different scan sizes.', 'err'); return;
    }
    msg('Processing …');
    setTimeout(() => {
      try {
        const t0 = performance.now();
        const { H, W } = fp, it = S.user.iters;
        const lp = level(fp.z, H, W, it), lm = level(fm.z, H, W, it);
        const A = gaussian(lp.mask.map((v) => 1 - v), H, W, 2), B = gaussian(lm.mask.map((v) => 1 - v), H, W, 2);
        const auto = crossCorrShift(A, B, H, W, 10);
        S.user.base = { fp, fm, lp, lm, auto, H, W, px: fp.pxNm / 1000, py: (fp.pyNm || fp.pxNm) / 1000, ms: performance.now() - t0 };
        composeUser();
        msg('');
      } catch (e) { msg('Processing failed: ' + e.message, 'err'); }
    }, 20);
  }

  function composeUser() {
    const b = S.user.base;
    if (!b) return;
    const { H, W, lp, lm } = b;
    const sh = S.user.manual || b.auto;
    const M2 = shiftBilinear(lm.z, H, W, sh[0], sh[1]);
    const D = new Float64Array(H * W);
    for (let i = 0; i < H * W; i++) D[i] = M2[i] - lp.z[i];
    const vr = validRegion(H, W, sh[0], sh[1]);
    if (vr.y1 - vr.y0 < 8 || vr.x1 - vr.x0 < 8) { msg('Shift too large: no overlap left.', 'err'); return; }
    const h0 = vr.y1 - vr.y0, w0 = vr.x1 - vr.x0;
    let P = crop(lp.z, H, W, vr), M = crop(M2, H, W, vr), Dc = crop(D, H, W, vr), patch = crop(lp.mask.map((v) => 1 - v), H, W, vr);
    const o = (a) => orient(a, h0, w0, S.user.flip, S.user.rot);
    const oP = o(P); P = oP.a; M = o(M).a; Dc = o(Dc).a; patch = o(patch).a;
    let h = oP.H, w = oP.W;
    let px = b.px, py = b.py;
    if (S.user.rot % 2) { const t = px; px = py; py = t; }
    if (S.user.crop) {
      const c = S.user.crop;
      P = crop(P, h, w, c); M = crop(M, h, w, c); Dc = crop(Dc, h, w, c); patch = crop(patch, h, w, c);
      h = c.y1 - c.y0; w = c.x1 - c.x0;
    }
    const mv = (f) => Math.round(Math.abs(f.bias) * 1000);
    const lab = (f, s) => (isFinite(f.bias) ? `${s}${mv(f)} mV` : `${s}V`);
    const lp2 = lab(b.fp, '+'), lm2 = lab(b.fm, '−');
    // step readouts
    const subPct = (m) => (100 * m.reduce((s, v) => s + v, 0) / m.length).toFixed(0);
    $('#cihdSteps').hidden = false;
    $('#cihdLevelOut').textContent = `substrate pixels: ${subPct(lp.mask)} % (+V), ${subPct(lm.mask)} % (−V) · ${b.ms.toFixed(0)} ms`;
    $('#cihdShiftAuto').textContent = `automatic shift: Δy = ${b.auto[0].toFixed(1)} px (${(b.auto[0] * b.py * 1000).toFixed(1)} nm), ` +
      `Δx = ${b.auto[1].toFixed(1)} px (${(b.auto[1] * b.px * 1000).toFixed(1)} nm)` + (S.user.manual ? ' · manual shift in use' : '');
    $('#cihdDy').value = sh[0].toFixed(1); $('#cihdDx').value = sh[1].toFixed(1);
    $('#cihdOverlap').textContent = `overlap ${w0} × ${h0} px (${(w0 * b.px).toFixed(2)} × ${(h0 * b.py).toFixed(2)} µm) of ${W} × ${H}`;
    $('#cihdRotOut').textContent = `${S.user.rot * 90}°${S.user.flip ? ', mirrored' : ''}`;
    $('#cihdCropOut').textContent = S.user.crop ? `crop ${w} × ${h} px (${(w * px).toFixed(2)} × ${(h * py).toFixed(2)} µm)` : 'no crop';
    const prev = S.view && S.view.source === 'user' ? S.view : null;
    const row = prev && prev.H === h ? prev.row : Math.floor(h / 2);
    setView({
      name: 'user', source: 'user', P, M, D: Dc, H: h, W: w, px, py, patch,
      labels: { plus: lp2, minus: lm2, plusCsv: lp2.replace('+', 'plus').replace(' ', ''), minusCsv: lm2.replace('−', 'minus').replace(' ', '') },
      row, half: prev ? prev.half : 2, preOrient: { h: h0, w: w0 },
      caption: `Your data: ${esc(b.fp.name)} (${lp2}) and ${esc(b.fm.name)} (${lm2}). Δh = h(${lm2}) − h(${lp2}) on the overlap, ` +
               `relative to the levelled substrate (= 0). Shift Δy = ${sh[0].toFixed(1)} px, Δx = ${sh[1].toFixed(1)} px.`,
    });
  }

  /* ---------------- shared viewer ---------------- */
  function setView(v) {
    const lim = (a, lo, hi) => { const s = sortedCopy(a); return [pctSorted(s, lo), pctSorted(s, hi)]; };
    v.limP = lim(v.P, 1, 99.5); v.limM = lim(v.M, 1, 99.5); v.limD = lim(v.D, 0.5, 99.5);
    v.imgP = mapImage(v.P, v.H, v.W, v.limP[0], v.limP[1], 'afmhot');
    v.imgM = mapImage(v.M, v.H, v.W, v.limM[0], v.limM[1], 'afmhot');
    v.imgD = mapImage(v.D, v.H, v.W, v.limD[0], v.limD[1], 'gray');
    v.statA = patchStatsP75(v.P, v.D, v.H, v.W);
    v.statB = patchStatsMask(v.patch, v.D, v.H, v.W);
    v.row = Math.max(0, Math.min(v.H - 1, v.row));
    S.view = v;
    $('#cihdView').hidden = false;
    const r = $('#cihdRow'); r.max = String(v.H - 1); r.value = String(v.row);
    $('#cihdHalf').value = String(v.half);
    renderMaps();
    renderProfile();
    $('#cihdCaption').innerHTML = v.caption;
    renderSummary();
  }

  function mapSpec(v, key, th) {
    const cb = { afmhot: cbarImage('afmhot'), gray: cbarImage('gray') };
    if (key === 'P') return { title: `Height at ${v.labels.plus}`, img: v.imgP, cbar: cb.afmhot, vmin: v.limP[0], vmax: v.limP[1], unit: 'Height (nm)', badge: v.labels.plus, badgeBg: COL.plus, badgeFg: '#fff', lineColor: COL.plus };
    if (key === 'M') return { title: `Height at ${v.labels.minus}`, img: v.imgM, cbar: cb.afmhot, vmin: v.limM[0], vmax: v.limM[1], unit: 'Height (nm)', badge: v.labels.minus, badgeBg: COL.minus, badgeFg: '#fff', lineColor: COL.minus };
    return { title: 'Δh = h(−V) − h(+V)', img: v.imgD, cbar: cb.gray, vmin: v.limD[0], vmax: v.limD[1], unit: 'Δh (nm)', badge: 'Δh', badgeBg: '#ffffff', badgeFg: '#111', lineColor: COL.dh };
  }

  function renderMaps() {
    const v = S.view, th = themeColors(false);
    root.querySelectorAll('.cihd-map').forEach((fig) => {
      const key = fig.dataset.key;
      const spec = Object.assign(mapSpec(v, key, th), { H: v.H, W: v.W, px: v.px, py: v.py, row: v.row, half: v.half });
      fig.innerHTML = mapSVG(spec, th);
      attachMapDrag(fig.querySelector('svg'));
    });
  }

  function svgPoint(svg, e) {
    const pt = svg.createSVGPoint(); pt.x = e.clientX; pt.y = e.clientY;
    return pt.matrixTransform(svg.getScreenCTM().inverse());
  }
  function attachMapDrag(svg) {
    const v = S.view, g = mapGeom(v.H, v.W, v.px, v.py);
    const toRC = (e) => {
      const p = svgPoint(svg, e);
      return { r: Math.floor(((p.y - g.MT) / g.PH) * S.view.H), c: Math.floor(((p.x - g.ML) / g.PW) * S.view.W) };
    };
    let drag = null;
    svg.addEventListener('pointerdown', (e) => {
      const rc = toRC(e);
      if (rc.r < 0 || rc.r >= S.view.H || rc.c < 0 || rc.c >= S.view.W) return;
      svg.setPointerCapture(e.pointerId);
      e.preventDefault();
      if (S.mode === 'user' && S.user.cropMode) { drag = { crop: true, r0: rc.r, c0: rc.c }; return; }
      drag = { crop: false };
      setRow(rc.r);
    });
    svg.addEventListener('pointermove', (e) => {
      if (!drag) return;
      const rc = toRC(e);
      if (drag.crop) { drag.r1 = rc.r; drag.c1 = rc.c; drawCropPreview(drag); return; }
      setRow(rc.r);
    });
    const end = (e) => {
      if (drag && drag.crop && isFinite(drag.r1)) {
        const v2 = S.view, cl = (x, n) => Math.max(0, Math.min(n, x));
        const y0 = cl(Math.min(drag.r0, drag.r1), v2.H), y1 = cl(Math.max(drag.r0, drag.r1) + 1, v2.H);
        const x0 = cl(Math.min(drag.c0, drag.c1), v2.W), x1 = cl(Math.max(drag.c0, drag.c1) + 1, v2.W);
        if (y1 - y0 >= 8 && x1 - x0 >= 8) {
          const base = S.user.crop || { y0: 0, x0: 0 };
          S.user.crop = { y0: base.y0 + y0, y1: base.y0 + y1, x0: base.x0 + x0, x1: base.x0 + x1 };
          S.user.cropMode = false;
          $('#cihdCropMode').setAttribute('aria-pressed', 'false');
          root.classList.remove('cihd-cropping');
          composeUser();
        } else { msg('Crop box too small (minimum 8 × 8 px).', 'err'); renderMaps(); }
      }
      drag = null;
    };
    svg.addEventListener('pointerup', end);
    svg.addEventListener('pointercancel', () => { drag = null; });
  }
  function drawCropPreview(d) {
    const v = S.view, g = mapGeom(v.H, v.W, v.px, v.py);
    root.querySelectorAll('.cihd-map svg').forEach((svg) => {
      let r = svg.querySelector('.cihd-croprect');
      if (!r) { r = document.createElementNS('http://www.w3.org/2000/svg', 'rect'); r.setAttribute('class', 'cihd-croprect'); r.setAttribute('fill', 'none'); r.setAttribute('stroke', '#ffffff'); r.setAttribute('stroke-width', '1.2'); r.setAttribute('stroke-dasharray', '3 2'); svg.appendChild(r); }
      const y0 = Math.min(d.r0, d.r1), y1 = Math.max(d.r0, d.r1) + 1, x0 = Math.min(d.c0, d.c1), x1 = Math.max(d.c0, d.c1) + 1;
      r.setAttribute('x', g.ML + (x0 / v.W) * g.PW); r.setAttribute('y', g.MT + (y0 / v.H) * g.PH);
      r.setAttribute('width', ((x1 - x0) / v.W) * g.PW); r.setAttribute('height', ((y1 - y0) / v.H) * g.PH);
    });
  }

  function setRow(r) {
    const v = S.view;
    if (!v) return;
    r = Math.max(0, Math.min(v.H - 1, r | 0));
    if (r === v.row && $('#cihdRow').value === String(r)) return;
    v.row = r;
    $('#cihdRow').value = String(r);
    updateLine();
  }
  function updateLine() {
    const v = S.view, g = mapGeom(v.H, v.W, v.px, v.py);
    const y = g.MT + ((v.row + 0.5) / v.H) * g.PH;
    const y0 = Math.max(0, v.row - v.half), y1 = Math.min(v.H - 1, v.row + v.half);
    root.querySelectorAll('.cihd-map svg').forEach((svg) => {
      const l = svg.querySelector('.cihd-line'), b = svg.querySelector('.cihd-band');
      l.setAttribute('y1', y); l.setAttribute('y2', y);
      b.setAttribute('y', g.MT + (y0 / v.H) * g.PH); b.setAttribute('height', ((y1 - y0 + 1) / v.H) * g.PH);
    });
    $('#cihdRowOut').textContent = `y = ${((v.row + 0.5) * v.py).toFixed(2)} µm (line ${v.row})`;
    renderProfile();
  }

  function currentProfile(v) {
    const x = Array.from({ length: v.W }, (_, i) => i * v.px);
    return {
      x, labels: v.labels,
      plus: Array.from(profileRows(v.P, v.H, v.W, v.row, v.half)),
      minus: Array.from(profileRows(v.M, v.H, v.W, v.row, v.half)),
      dh: Array.from(profileRows(v.D, v.H, v.W, v.row, v.half)),
    };
  }

  function renderProfile() {
    const v = S.view;
    if (!v) return;
    const box = $('#cihdProf');
    const p = currentProfile(v);
    const width = box.clientWidth || 900;
    const narrow = width < 520;
    const { svg, geom } = profileSVG(p, narrow ? width * 1.25 : width, themeColors(false), { height: narrow ? 300 : 270 });
    box.innerHTML = svg;
    $('#cihdRowOut').textContent = `y = ${((v.row + 0.5) * v.py).toFixed(2)} µm (line ${v.row})`;
    const el = box.querySelector('svg'), hit = el.querySelector('.cihd-hit'), hov = el.querySelector('.cihd-hover');
    const ln = hov.querySelector('line'), tx = hov.querySelector('text');
    hit.addEventListener('pointermove', (e) => {
      const q = svgPoint(el, e);
      const i = Math.max(0, Math.min(p.x.length - 1, Math.round(((q.x - geom.ML) / geom.PW) * geom.x1 / v.px)));
      const X = geom.ML + (p.x[i] / geom.x1) * geom.PW;
      ln.setAttribute('x1', X); ln.setAttribute('x2', X);
      const left = X > geom.ML + geom.PW * 0.62;
      tx.setAttribute('x', left ? X - 6 : X + 6); tx.setAttribute('y', 44);
      tx.setAttribute('text-anchor', left ? 'end' : 'start');
      tx.textContent = `${p.x[i].toFixed(3)} µm   ${p.labels.plus} ${p.plus[i].toFixed(2)}   ${p.labels.minus} ${p.minus[i].toFixed(2)}   Δh ${p.dh[i].toFixed(2)} nm`;
      hov.style.display = '';
    });
    hit.addEventListener('pointerleave', () => { hov.style.display = 'none'; });
  }

  function renderSummary() {
    const v = S.view, A = v.statA, B = v.statB;
    const rowTxt = (name, s, note) => `<tr><td>${name}</td><td class="num">${signed(s.median, 1)}</td><td class="num">${signed(s.q25, 1)} to ${signed(s.q75, 1)}</td><td class="num">${s.n}</td><td>${note}</td></tr>`;
    let html = `<table class="cihd-table"><caption>Δh on the patches (nm), relative to the substrate</caption>
      <thead><tr><th>Patch definition</th><th>Median</th><th>IQR</th><th>Pixels</th><th></th></tr></thead><tbody>`;
    html += rowTxt('Levelling mask, 2× eroded', B, 'default');
    html += rowTxt('+V height > 75th percentile, 2× eroded', A, 'as in the Python reference');
    if (v.published) {
      const p = v.published;
      html += `<tr><td>Published (re-processed)</td><td class="num">${signed(p.median, 1)}</td><td class="num">${signed(p.q25, 1)} to ${signed(p.q75, 1)}</td><td class="num">—</td><td>README</td></tr>`;
    }
    html += '</tbody></table>';
    const sgn = B.median > 0.5 ? 'positive' : B.median < -0.5 ? 'negative' : 'close to zero';
    html += `<p class="cihd-note">Δh on the patches is ${sgn}: the patches are ${B.median > 0.5 ? 'more positive' : B.median < -0.5 ? 'more negative' : 'not clearly different in charge'} than the substrate. ` +
            `The 75th-percentile rule depends on the field of view (crop), so the two values can differ.</p>`;
    $('#cihdSummary').innerHTML = html;
  }

  /* ---------------- export ---------------- */
  function download(name, data, type) {
    const blob = data instanceof Blob ? data : new Blob([data], { type });
    const a = document.createElement('a');
    a.href = URL.createObjectURL(blob); a.download = name;
    document.body.appendChild(a); a.click(); a.remove();
    setTimeout(() => URL.revokeObjectURL(a.href), 2000);
  }
  function matrixCSV(a, H, W) {
    const lines = [];
    for (let y = 0; y < H; y++) { const r = []; for (let x = 0; x < W; x++) r.push(a[y * W + x].toFixed(3)); lines.push(r.join(',')); }
    return lines.join('\n') + '\n';
  }
  function figureSVG() {
    const v = S.view, th = themeColors(true);
    const FW = 1200, gap = 16, mw = (FW - 2 * gap - 40) / 3;
    const g = mapGeom(v.H, v.W, v.px, v.py), sc = mw / g.VW, mh = g.VH * sc;
    let s = `<svg xmlns="http://www.w3.org/2000/svg" width="${FW}" height="__H__" viewBox="0 0 ${FW} __H__" font-family="Geist, Helvetica, Arial, sans-serif"><rect width="100%" height="100%" fill="#ffffff"/>`;
    ['P', 'M', 'D'].forEach((k, j) => {
      const spec = Object.assign(mapSpec(v, k, th), { H: v.H, W: v.W, px: v.px, py: v.py, row: v.row, half: v.half });
      const inner = mapSVG(spec, th).replace('<svg ', `<svg x="${(20 + j * (mw + gap)).toFixed(1)}" y="10" width="${mw.toFixed(1)}" height="${mh.toFixed(1)}" `);
      s += inner;
    });
    const pr = profileSVG(currentProfile(v), FW - 40, th, { height: 300 });
    const py0 = 10 + mh + 10;
    s += pr.svg.replace('<svg ', `<svg x="20" y="${py0.toFixed(1)}" width="${FW - 40}" height="300" `);
    const cap = (v.caption || '').replace(/<[^>]+>/g, '');
    const words = cap.split(' '); let line = '', ly = py0 + 300 + 22; const lines = [];
    for (const w of words) { if ((line + ' ' + w).length > 180) { lines.push(line); line = w; } else line = line ? line + ' ' + w : w; }
    if (line) lines.push(line);
    lines.forEach((l) => { s += `<text x="20" y="${ly}" font-size="12" fill="#333">${esc(l)}</text>`; ly += 16; });
    const H = Math.ceil(ly + 6);
    return s.replace(/__H__/g, String(H)) + '</svg>';
  }
  function doExport(kind) {
    const v = S.view;
    if (!v) return;
    const base = v.source === 'example' ? `CIHD_${v.name}` : 'CIHD_user';
    if (kind === 'maps') {
      download(`${base}_map_${v.labels.plusCsv}_nm.csv`, matrixCSV(v.P, v.H, v.W), 'text/csv');
      setTimeout(() => download(`${base}_map_${v.labels.minusCsv}_nm.csv`, matrixCSV(v.M, v.H, v.W), 'text/csv'), 300);
      setTimeout(() => download(`${base}_map_dh_nm.csv`, matrixCSV(v.D, v.H, v.W), 'text/csv'), 600);
    } else if (kind === 'profile') {
      const p = currentProfile(v);
      let t = `position_um,height_${v.labels.plusCsv}_nm,height_${v.labels.minusCsv}_nm,dh_nm (dh = h(${v.labels.minus.replace('−', '-')}) - h(${v.labels.plus})); mean of ${2 * v.half + 1} lines around y = ${((v.row + 0.5) * v.py).toFixed(4)} um\n`;
      for (let i = 0; i < p.x.length; i++) t += `${p.x[i].toFixed(4)},${p.plus[i].toFixed(4)},${p.minus[i].toFixed(4)},${p.dh[i].toFixed(4)}\n`;
      download(`${base}_profile.csv`, t, 'text/csv');
    } else if (kind === 'svg') {
      download(`${base}_figure.svg`, figureSVG(), 'image/svg+xml');
    } else if (kind === 'png') {
      const svg = figureSVG();
      const W = +svg.match(/width="(\d+)"/)[1], H = +svg.match(/height="(\d+)"/)[1];
      const img = new Image();
      img.onload = () => {
        const cv = document.createElement('canvas'); cv.width = 2 * W; cv.height = 2 * H;
        const ctx = cv.getContext('2d'); ctx.scale(2, 2); ctx.drawImage(img, 0, 0);
        cv.toBlob((b) => download(`${base}_figure.png`, b), 'image/png');
      };
      img.onerror = () => msg('PNG export failed in this browser; use SVG.', 'err');
      img.src = 'data:image/svg+xml;charset=utf-8,' + encodeURIComponent(svg);
    }
  }

  /* ---------------- init ---------------- */
  function init() {
    root = document.getElementById('cihd-root');
    if (!root) return;
    buildUI();
    // example mode loads when the section is first shown (keeps the landing page light)
    let started = false;
    const start = () => { if (!started) { started = true; setMode('example'); } };
    const tab = document.getElementById('content-qscm');
    if (tab && window.MutationObserver) {
      const ob = new MutationObserver(() => { if (tab.style.display !== 'none' && tab.classList.contains('active')) start(); });
      ob.observe(tab, { attributes: true, attributeFilter: ['style', 'class'] });
      if (tab.classList.contains('active') && tab.style.display !== 'none') start();
    } else start();
    if (location.hash === '#cihd' && typeof window.openSuite === 'function') {
      setTimeout(() => { window.openSuite('qscm'); start(); root.scrollIntoView({ block: 'start' }); }, 60);
    }
  }
  if (document.readyState === 'loading') document.addEventListener('DOMContentLoaded', init); else init();

  window.CIHD = {
    core: { percentile, median, stdPop, level, gaussian, crossCorrShift, shiftBilinear, validRegion, crop, orient, erode,
            patchStatsP75, patchStatsMask, maskedStats, profileRows, pyRound },
    readParkImage, readTextMatrix, loadUserFiles, setMode, loadExample, setRow,
    state: S,
    // for tests: export contents without downloading
    exportContent: { figureSVG, matrixCSV: () => S.view && matrixCSV(S.view.D, S.view.H, S.view.W), profile: () => S.view && currentProfile(S.view) },
  };
})();
