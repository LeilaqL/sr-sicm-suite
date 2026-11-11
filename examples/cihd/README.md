# CIHD example – lipid patches on mica (dual-bias SICM imaging)

Data: six Park XEP topography images (2015), DPTAP (cationic), DPPE (zwitterionic), DPPG (anionic)
patches on mica, each imaged at +100 mV and −100 mV. Source: Klausen et al., Nat. Commun. 7, 12447 (2016)
(dataset shared on request – cite the paper on the website).

Processing (py/ folder):
1. Height = raw × data gain (Park header tag 50435) → nm.
2. Levelling on the mica substrate: plane fit + line-by-line median, iterating to exclude lipid patches; substrate set to 0.
3. Registration of the −100 mV image onto the +100 mV image by cross-correlating the patch masks (sub-pixel), shift applied.
4. Δh = h(−100 mV) − h(+100 mV)  (same sign rule as Δz99: Δh > 0 → more positive than the mica substrate).
5. Orientation and crop as in the published Fig. 2: the images were mirrored/rotated (DPTAP: flip + 180°, DPPE and DPPG: flip + 270°)
   and cropped by template matching against the published −100 mV panels: DPTAP and DPPE 4 × 2 µm, DPPG 1.8 × 0.9 µm.
   Profiles: mean of 5 scan lines at the dashed line (y = 1.10, 0.92 and 0.40 µm, as in the paper).

Bias label: the paper's "+100 mV" corresponds to header "sample bias" −0.10 V in these 2015 files
(pipette potential = −header value). Labels in the figures follow the paper.

Result (median Δh on patches, IQR): DPTAP +9.2 nm (6.2–11.3), DPPE +7.0 nm (6.5–7.4), DPPG −2.1 nm (−3.4 to −1.2).
Δh is relative to the (negative) mica substrate.

Caption (EN): "Dual-bias SICM images of lipid patches on mica. Topography at +100 mV (top) and −100 mV (middle),
their difference Δh = h(−100 mV) − h(+100 mV) (third row) and line profiles along the dashed lines (bottom).
Cationic DPTAP and zwitterionic DPPE appear higher at −100 mV (Δh > 0), anionic DPPG lower (Δh < 0):
the sign of Δh follows the surface charge relative to the substrate. Data: Klausen et al., Nat. Commun. 2016; re-processed."

Files: CIHD_lipids_all.(png|svg), CIHD_<lipid>.(png|svg), CIHD_summary_bar.(png|svg),
CIHD_<lipid>_profile.csv, CIHD_<lipid>_map_{plus100mV,minus100mV,dh}_nm.csv (registered, oriented, cropped, nm; pixel size 19.5 nm for DPTAP/DPPE, 7.8 nm for DPPG).
