# Simulated demo: dual-bias triplet readout

**Simulated data, not measured.** Used for the "Load example" button of the triplet module and as a validation test:
the site must recover the known Δz99.

## What is simulated
- Δz99 = z99(V−) − z99(V+) = **−6.0 nm** (negative surface).
- Surface drift −14 nm/min with slight curvature (+0.04 nm/min²).
- 10 blocks alternating +0.1 V / −0.1 V, 5 curves per block, 14 s apart, one block every 3 min.
- The first 2 curves of each block are displaced by +12 and +5 nm (settling after the bias switch); the analysis discards them.
- Block 7: I∞ changes by 1 % over its 5 curves (0.5 % over the 3 used curves), so it must be flagged (criterion 0.3 %).
- Current noise 0.03 % of I∞. Curve shape I/I∞ = 1 − 0.01 exp(−(h − z99)/20 nm).

## Reference result (Python, `sicm_lib`, unchanged)
z99 from the bulk side at 0.99, blocks with the first 2 curves discarded, triplets with times in seconds:

| | n | Δz99 mean ± SD (nm) |
|---|---|---|
| all triplets | 8 | −6.19 ± 0.72 |
| without triplets that use the flagged block 7 | 5 | −6.31 ± 0.78 |
| for comparison: −V block minus the previous +V block, no drift correction | 5 | −44.8 |

Full numbers in `expected.json`. Regenerate with `python gen_demo.py <path to sicm_tools>` (fixed random seed).

## Files
- `demo_curves.json`: curves (file name, block, repetition, time, bias, start height, step, current in nA).
- `expected.json`: blocks, triplets and summary from the Python reference.
- `gen_demo.py`: generator.
