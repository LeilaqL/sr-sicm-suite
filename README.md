# SR-SICM Suite

Browser-based analysis tools for scanning ion conductance microscopy (SICM): super-resolution imaging and surface-charge (QSCM) analysis of Park Systems XE/XEI data.

**Live site:** https://sr-sicm-suite.netlify.app


## What it does

| Module | Input | Output |
|---|---|---|
| Surface charge – dual-bias approach curves | Park TIFF approach curves at ±V | Δz₉₉ per drift-corrected triplet, mean ± SD |
| Surface charge – dual-bias imaging (CIHD) | Two topography images (+V, −V) of the same area | Levelled and registered maps, Δh map, interactive line profiles |
| Asymmetric conductivity | I–V curves | Rectification and offset V₀ |
| Z-detector calibration | Calibration grid scan | nm per count |
| Stiffness | Approach curves on sample and rigid reference | Slope ratio s∞/s, indentation δ |

All processing runs locally in the browser; no data are uploaded to a server.

## Conventions

- Δz₉₉ = z₉₉(V−) − z₉₉(V+); its sign is the sign of the surface charge.
- z₉₉: position where I/I∞ = 0.99, approached from the bulk side.
- Slope s: linear fit between 0.98 and 0.99 I/I∞, per curve.
- Blocks of 5 curves per bias, first 2 discarded; I∞ drift > 0.3 % within a block is flagged.
- Z detector: 0.879 nm/count (calibrated on a 1.5 µm PDMS grid).

## Example data

`examples/cihd/` contains processed maps of DPTAP, DPPE and DPPG lipid patches on mica (±100 mV), re-processed from Klausen et al., *Nat. Commun.* (2016), shown with permission of the data owner. Raw data are not included.

## Run locally

```
python -m http.server 8765
```
then open http://127.0.0.1:8765 (the examples do not load from `file://`).

## References

- Rheinlaender J., Schäffer T.E., *Soft Matter* 9, 3230 (2013) – stiffness from approach-curve slopes.
- Klausen L.H., Fuhs T., Dong M., *Nat. Commun.* (2016) – dual-bias surface-charge imaging.
- Skaanvik et al. (2025) – asymmetric conductivity.

## Author

Leyla Nur Akay, Department of Chemistry, Aarhus University (supervisor: Prof. Mingdong Dong).

## License

Code: MIT. Example data: shown with permission; not for reuse without the data owner's consent.
