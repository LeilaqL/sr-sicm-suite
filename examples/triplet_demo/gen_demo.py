"""Simulated demo data for the dual-bias triplet readout (SR-SICM Suite).

NOT MEASURED DATA. Approach curves are generated with a known answer:
  * dz99 = z99(V-) - z99(V+) = -6.0 nm (negative surface)
  * surface drift z_s(t) = -14 nm/min * t + 0.04 nm/min^2 * t^2 (slightly curved)
  * 10 blocks alternating +0.1 V / -0.1 V, 5 curves per block, 14 s apart, one block every 3 min
  * the first 2 curves of each block are displaced by +12 and +5 nm (settling after the bias switch)
  * block 7 has an I_inf drift of 1 % over its 5 curves, 0.5 % over the 3 used curves (must be flagged, criterion 0.3 %)
  * current noise 0.03 % of I_inf
Curve shape: I/I_inf = 1 - 0.01 exp(-(h - z99)/L), L = 20 nm, so I/I_inf = 0.99 exactly at h = z99.

Writes demo_curves.json (input for the site) and expected.json (result of the unchanged Python
reference sicm_lib: z_cross, blocks(discard=2), triplets).
Usage: python gen_demo.py <path to sicm_tools>
"""
import sys, json, datetime as dt
import numpy as np

sys.path.insert(0, sys.argv[1])
import sicm_lib as L  # reference implementation, not modified

TRUE_DZ99 = -6.0
DRIFT_LIN, DRIFT_QUAD = -14.0, 0.04          # nm/min, nm/min^2
NBLOCK, NCURVE, REP_S, BLOCK_S = 10, 5, 14, 180
LEN_NM, STEP_NM, DECAY_NM = 1200.0, 3.0, 20.0
SETTLE = [12.0, 5.0, 0.0, 0.0, 0.0]          # nm, first two curves of each block
NOISE = 3e-4
T0 = dt.datetime(2026, 1, 1, 10, 0, 0)
rng = np.random.default_rng(20261111)

def z_surface(t_min):
    return 5000.0 + DRIFT_LIN * t_min + DRIFT_QUAD * t_min ** 2

curves_out, ref_curves = [], []
for b in range(NBLOCK):
    V = 0.1 if b % 2 == 0 else -0.1
    I_inf0 = (2.000 if V > 0 else 1.800) * (1 + 0.0005 * b)
    for r in range(NCURVE):
        t = T0 + dt.timedelta(seconds=b * BLOCK_S + r * REP_S)
        tm = (t - T0).total_seconds() / 60
        z99 = z_surface(tm) + (TRUE_DZ99 / 2 if V < 0 else -TRUE_DZ99 / 2) + SETTLE[r]
        I_inf = I_inf0 * (1 + (0.010 * r / (NCURVE - 1) if b == 6 else 0.0))   # 0.5 % over the 3 used curves
        h = z99 + LEN_NM - STEP_NM * np.arange(int((LEN_NM + 30) / STEP_NM))   # approach: far -> near
        rel = 1 - 0.01 * np.exp(-(h - z99) / DECAY_NM)
        I = np.sign(V) * I_inf * (rel + NOISE * rng.standard_normal(h.size))
        name = f'demo_block{b + 1:02d}_{r + 1}.tiff'
        curves_out.append(dict(file=name, block=b + 1, rep=r + 1, time=t.isoformat(), bias_V=V,
                               h0_nm=round(float(h[0]), 3), dh_nm=-STEP_NM, I_nA=[round(float(x), 6) for x in I]))
        # the same curve as sicm_lib.load_curves would build it (h ascending, I normalised by the far median)
        far = (h - h[-1]) > 1000
        I0 = float(np.median(np.abs(I)[far]))
        s = np.argsort(h)
        c = dict(file=name, rep=1, h=h[s], I=np.abs(I)[s] / I0, I0=I0, bias=V, time=t, zturn=h[-1])
        ref_curves.append(L.analyse(c))

B = L.blocks(ref_curves, discard=2)
T = L.triplets(B)
d = np.array([t['dz99'] for t in T])
flag = [i + 1 for i, b in enumerate(B) if b['I0_drift_pct'] > 0.3]
keep = [t for t in T if t['max_I0_drift_pct'] <= 0.3]
dk = np.array([t['dz99'] for t in keep])
# naive comparison without drift correction: each -V block against the previous +V block
naive = [B[k]['z99'] - B[k - 1]['z99'] for k in range(1, len(B)) if B[k]['bias'] < 0]

expected = dict(
    note='Simulated demo data, not measured. Reference result of sicm_lib (unchanged): z_cross at 0.99 from the bulk side, '
         'blocks(discard=2), triplets(). Times from the curve time stamps (seconds).',
    true_dz99_nm=TRUE_DZ99, drift_nm_per_min=DRIFT_LIN, drift_curvature_nm_per_min2=DRIFT_QUAD,
    blocks=[dict(index=i + 1, bias_V=b['bias'], n_used=b['n'], z99_nm=round(b['z99'], 3),
                 I_inf_nA=round(b['I0'], 5), I_inf_drift_pct=round(b['I0_drift_pct'], 3)) for i, b in enumerate(B)],
    triplets=[dict(blocks=[t['blocks'][0], t['blocks'][1], t['blocks'][2]], dz99_nm=round(t['dz99'], 3),
                   max_I_inf_drift_pct=round(t['max_I0_drift_pct'], 3)) for t in T],
    dz99_all=dict(n=len(d), mean_nm=round(float(d.mean()), 3), sd_nm=round(float(d.std(ddof=1)), 3)),
    dz99_without_flagged=dict(n=len(dk), mean_nm=round(float(dk.mean()), 3), sd_nm=round(float(dk.std(ddof=1)), 3)),
    flagged_blocks=flag,
    naive_pairs_without_drift_correction=dict(n=len(naive), mean_nm=round(float(np.mean(naive)), 3)),
)
json.dump(dict(title='Simulated demo: dual-bias triplet series (not measured data)',
               description='10 blocks, +/-0.1 V, 5 curves per block (14 s apart), known dz99 = -6.0 nm, '
                           'drift -14 nm/min with slight curvature, current noise 0.03 %, block 7 with 0.5 % I_inf drift over the used curves.',
               simulated=True, curves=curves_out), open('demo_curves.json', 'w'), separators=(',', ':'))
json.dump(expected, open('expected.json', 'w'), indent=1)
print(json.dumps({k: expected[k] for k in ('dz99_all', 'dz99_without_flagged', 'flagged_blocks', 'naive_pairs_without_drift_correction')}, indent=1))
for t in expected['triplets']:
    print(t['blocks'][1], t['dz99_nm'], t['max_I_inf_drift_pct'])
