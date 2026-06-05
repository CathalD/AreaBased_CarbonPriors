# Area-Based Carbon Priors — Forest & Soil Carbon Prior Assessment and Field-Sampling Design

**Version:** v4.6 · **Prepared:** 2026 · **Pipeline:** `GEE_Script.js` (Google Earth Engine)
**Area of interest (AOI):** *set per project — any GEE FeatureCollection boundary loaded in Step 0 (carbon priors valid for AOIs within Canada)*
**Projection / scale:** UTM, auto-selected from the AOI centroid (`AUTO_UTM`; example run EPSG:32621 / UTM Zone 21N) · forest 25 m, soil 250 m native

> **Note on numbers.** This is a reusable report template. The figures below come
> from one example end-to-end run of `GEE_Script.js` v4.6 (`[▶ RUN ALL]`,
> `N_UNC_BINS = 3`, forest n = 10, soil n = 12) over a Canadian AOI; re-run the
> pipeline on your own boundary and substitute the printed values. AOI-mean stats
> use Earth Engine `bestEffort` reductions and may vary ~±0.2 kg/m² between runs
> (see §8); exports are deterministic.

---

## 1. Purpose

This work compiles the best available **prior** maps of forest and soil carbon over
the AOI, quantifies where they agree and disagree (including their stated
uncertainties), combines them into a single uncertainty-weighted best estimate, and
uses that — together with land cover — to design a **statistically optimal field
campaign**. The deliverables are:

1. Harmonised forest and soil carbon prior maps (kg C/m²).
2. Per-product uncertainty maps and an inter-product disagreement map.
3. A combined, inverse-variance-weighted carbon map with propagated uncertainty.
4. A Neyman-allocated stratified sampling design.
5. Exactly **`N_FOREST_SAMPLES`** forest and **`N_SOIL_SAMPLES`** soil sample
   locations (default 10 and 12) distributed as points.
6. Report-style exports (rasters, a summary CSV, and KML/CSV point files) the
   field partner can load directly into Google Earth or GIS.

The script is built to be **handed to a partner** to run themselves: a one-click
`[▶ RUN ALL]` button executes the whole pipeline, and a `[0] Check Asset Access`
button verifies up front that every required dataset is readable.

---

## 2. Study area

The project AOI, defined by the boundary asset loaded in Step 0
(`AOI_ASSET`, e.g. `projects/your-gee-project/assets/your_aoi_boundary`). The pipeline
buffers the AOI by 50 km (`AOI_BUFFER_M`) for sampling remote-sensing predictors, then
clips all outputs to the boundary.

- AOI area: **8,438.8 ha**
- Dominant land cover (ESA WorldCover): predominantly **forest**, with minor
  shrub/grassland; negligible mapped wetland/cropland within the AOI (the Neyman
  strata that received area were almost entirely the forest and shrub/grass classes).

---

## 3. Data sources

### 3.1 Forest carbon members (3)

| Member | Source | Native scope | Resolution | Uncertainty (σ) used |
|---|---|---|---|---|
| **GEDI RF** (this study) | Random-forest model trained on GEDI L4A AGBD, extrapolated with optical/SAR/terrain covariates | Aboveground biomass density (AGBD) | 25 m | 3-fold CV RMSE (constant) |
| **Sothe FC** | Sothe et al. (2022), McMaster / WWF-Canada national forest-carbon map | **Total vegetation carbon: aboveground + belowground (roots) + dead plants** | 250 m | Per-pixel published uncertainty layer |
| **SCANFI v1.2** | Canada SCANFI biomass (NRCan) | Aboveground biomass (AGB) | ~30 m | Constant, from Guindon (2024) RMSE 38.70 t/ha |

### 3.2 Soil carbon members (2)

| Member | Source | Scope | Resolution | Uncertainty (σ) used |
|---|---|---|---|---|
| **Sothe SC** | Sothe et al. (2022) national soil-carbon map | Soil organic carbon, 0–1 m | 250 m | Per-pixel published uncertainty layer |
| **SoilGrids OCS** | ISRIC SoilGrids v2.0 (SOC × bulk density, depth-integrated 0–100 cm) | Soil organic carbon, 0–1 m | 250 m | **Derived** (see §4.4) — no per-pixel layer published in GEE |

### 3.3 Covariates (forest RF)

NASADEM elevation/slope/aspect, TWI, TPI; Sentinel-2 summer-median bands + NDVI,
NDWI, NBR, EVI; Sentinel-1 SAR seasonal composites; TerraClimate (tmin, tmax,
precip, soil moisture); Canada SBFI structure AGB mean/SD; lead tree species.
Variable selection keeps the top `TOP_N_BANDS` (default 6) by Gini importance.

### 3.4 Optional field data (Step 1 / Step 4b only)

WOSIS 2023, CanPeat peat profiles, and combined soil profiles — used only for the
covariate-extraction CSVs, **not** for the core carbon/sampling pipeline.

---

## 4. Methods

### 4.1 Carbon harmonisation (the key v4.6 fix)

The three forest members are **not natively the same quantity**. Sothe FC is total
vegetation carbon including **below-ground (root) biomass**, whereas GEDI L4A and
SCANFI are **above-ground biomass only**. Averaging them as-is would mix biomass
with carbon and above-ground-only with above+below-ground.

All forest members are therefore converted to a common quantity — **above- + below-
ground biomass carbon (kg C/m²)** — using:

```
Carbon (kg C/m²) = AGB (Mg/ha)
                 × (1 + BGB_RATIO)     // add belowground: BGB = 0.22 × AGB
                 × CARBON_FRACTION     // dry-biomass carbon fraction = 0.50
                 × MGHA_TO_KGM2        // Mg/ha → kg/m² = 0.10
                 = AGB × 0.061
```

- `BGB_RATIO = 0.22` — belowground:aboveground biomass ratio (BGB ≈ 0.22 × AGB).
- `CARBON_FRACTION = 0.50` — IPCC default carbon fraction of dry biomass.
- `AGB_TO_C_KGM2 = 0.061` — combined factor applied to GEDI and SCANFI **and to
  their σ** (so weights stay in carbon units).

> **Residual scope difference (declared, not corrected).** Sothe additionally
> includes **dead-plant carbon** and used **forest-type-specific** root ratios,
> while our GEDI/SCANFI conversion uses a single national BGB ratio of 0.22. Sothe
> FC will therefore read slightly higher than the GEDI/SCANFI members where dead
> material is abundant. This is documented in the summary table's `carbon_scope`
> column and is small relative to the per-product uncertainties.

Soil members are already soil organic carbon over 0–1 m and need no biomass
conversion.

### 4.2 Forest carbon model (GEDI RF)

A two-stage random forest:

1. **Pilot RF** (`PILOT_N_TREES` trees) on all candidate covariates → Gini variable
   importance → keep top `TOP_N_BANDS`.
2. **Final RF** (`FINAL_N_TREES` trees) on the selected bands, trained on GEDI L4A
   AGBD footprints (quality- and degrade-filtered, `0 < AGBD < 600`), then applied
   wall-to-wall and clipped to plausible range.

**3-fold cross-validation** holds out one third of footprints at a time, trains on
the rest, and computes hold-out RMSE; the mean of the three folds is the GEDI RF
uncertainty (converted to carbon σ with `AGB_TO_C_KGM2`). This σ measures how well
the RF reproduces held-out GEDI values — it is *model* error given GEDI as truth,
not the absolute error in the true carbon stock.

- GEDI L4A training footprints: 110 (109 complete cases); AGBD mean 35.4 Mg/ha,
  SD 26.1 Mg/ha.
- 3-fold CV RMSE: **22.26 Mg/ha AGBD** (folds 27.70 / 18.04 / 21.02) → σ =
  **1.358 kg C/m²**.
- Top 6 covariates retained: **B3, NBR, B12, B11, B2, SAR ascending-summer VH**
  (Sentinel-2 green/SWIR/NIR and SAR dominate).

### 4.3 Soil carbon

Two independent products (Sothe SC, SoilGrids OCS 0–100 cm) are compared and
differenced; their pixel-wise disagreement (SD) feeds both the SoilGrids σ
derivation and the sampling stratification.

### 4.4 Inverse-variance-weighted ensembles

Each pool is combined per pixel by precision weighting:

```
w_i(p) = 1 / σ_i(p)²
wmean(p) = Σ[ w_i(p) · x_i(p) ] / Σ[ w_i(p) ]
σ_wmean(p) = 1 / sqrt( Σ[ w_i(p) ] )
```

- **Forest (3 members)** — GEDI RF (constant σ), Sothe FC (per-pixel σ), SCANFI
  (constant σ). v4.6: a member that is missing at a pixel contributes **zero
  weight** rather than nulling the whole pixel; the result is kept only where
  **≥2 of 3** members are present.
- **Soil (2 members)** — Sothe SC (per-pixel σ) and SoilGrids.
  - **Derived SoilGrids σ.** No per-pixel SoilGrids uncertainty layer is available
    in GEE, so it is estimated as: standardise the Sothe σ (mean = 1) and the
    inter-product disagreement SD (mean = 1), average them equally, then rescale to
    the Sothe σ mean. Rationale: SoilGrids is most uncertain where Sothe is
    uncertain **and** where the two products disagree; standardising first prevents
    either component dominating by magnitude.

An **equal-weight** mean (same members) is retained for comparison. v4.6 aligned the
equal-weight forest ensemble to the same 3 members as the weighted one, so
"weighted vs equal-weight" is a like-for-like comparison.

### 4.5 Map agreement / disagreement

The pipeline produces difference layers (Sothe − SoilGrids for soil; Sothe − SCANFI
and weighted − equal-weight for forest), a forest data-source map (number of
products present, 1–3), and a soil data-source map. Summary agreement:

- Forest Sothe − SCANFI mean difference: **+1.51 kg C/m²** (range −7.9 to +13.3);
  Sothe FC reads higher than SCANFI on average.
- Soil Sothe − SoilGrids mean difference: **+9.17 kg C/m²** (range −36.7 to +84.4);
  Sothe SC reads substantially higher than SoilGrids, with large local disagreement.
- Weighted − equal-weight: forest **−0.98 kg C/m²**, soil **+0.28 kg C/m²**. The
  forest weighted mean sits well below the equal-weight mean because GEDI RF (lowest
  σ → highest weight) pulls the estimate down (see §8, GEDI weighting caveat).

### 4.6 Sampling design (Neyman allocation)

Strata = **land-cover class** × **uncertainty tertile bin**. Only two land-cover
classes are sampled — **Forest** (ESA WorldCover tree cover) and **Wetland**
(GWL_FCS30 swamp / marsh / flooded-flat / saline / mangrove / salt-marsh / tidal-flat;
open water and non-wetland excluded). The GWL_FCS30 wetland layer **overrides** the
ESA forest layer where they overlap, so treed peatlands are counted as wetland.
Shrub/grass and other cover are mapped for context but **not sampled**. The
uncertainty bins use `N_UNC_BINS = 3` (breaks at the 33rd/67th percentiles of the RSS
uncertainty); the binning generalises to any `N_UNC_BINS`.
Neyman optimal allocation then splits points across the Forest and Wetland strata in
proportion to **area × uncertainty** (`N_h × σ_h`), so high-uncertainty wetland — even
where it covers little area — draws samples.

v4.6 uses **largest-remainder rounding** so the per-stratum counts sum to **exactly**
`N_FOREST_SAMPLES` / `N_SOIL_SAMPLES` — whatever is set in the config is exactly what
is produced. If a stratum receives 0 or fewer than `MIN_PTS_PER_STRATUM` points, the
Console flags it (too few samples for the number of strata → raise N or reduce
`N_UNC_BINS`).

A **power analysis** reports, for each pool, the spatial CV, the SRS-bound minimum n
for a ±20% margin of error at 90% confidence, the expected MOE at the configured n,
and an explicit recommendation. It is an SRS bound (uses map spatial SD as a proxy
for plot variance), so Neyman stratification should do at least this well.

- Forest: μ = 3.386 kg/m², σ = 0.834, **CV = 24.6%**, **n_min(±20%) = 5**,
  configured n = 10 → expected MOE **±12.8%** → ✓ **comfortably meets the target**.
- Soil: μ = 38.277 kg/m², σ = 16.414, **CV = 42.9%**, **n_min(±20%) = 13**,
  configured n = 12 → expected MOE **±20.4%** → ⚠ **just misses** the ±20% SRS bound.
  **Recommendation: set `N_SOIL_SAMPLES = 13`** (Neyman stratification typically
  beats the SRS bound, so 12 may suffice in practice, but 13 guarantees ±20%).

---

## 5. Results

### 5.1 Carbon densities and stocks (weighted ensemble)

AOI area = **8,438.8 ha**.

| Pool | Mean density (Mg C/ha) | Mean density (kg C/m²) | Total stock (t C) |
|---|---|---|---|
| Forest (forested pixels) | 35.8 | 3.575 | **301,730** |
| Soil (0–1 m) | 382.8 | 38.277 | **3,230,090** |
| **Total ecosystem** | **418.5** | **41.852** | **3,531,820** |

> Total stock = mean density × AOI area. Soil holds ~91% of the ecosystem carbon
> prior (~3.23 Mt C of ~3.53 Mt C), consistent with Sothe et al.'s national finding
> that most Canadian ecosystem carbon is below ground. This is a **prior** from
> remote sensing; field sampling will produce the calibrated stock with a defensible
> CI.

### 5.2 Per-product comparison (within AOI)

| Source | Pool | Scope | σ type | σ mean (kg/m²) | Mean (kg/m²) | SD (kg/m²) |
|---|---|---|---|---|---|---|
| GEDI RF (this study) | Forest | AGB+BGB carbon | 3-fold CV | 1.358 | 2.83–3.09 † | 0.83–1.14 |
| Sothe et al. FC | Forest | AGB+BGB+dead | per-pixel | 4.721 | 5.981 | 2.345 |
| SCANFI v1.2 FC | Forest | AGB+BGB carbon | constant | 2.361 | 4.406 | 1.489 |
| Forest equal-weight mean | Forest | AGB+BGB carbon | ensemble SD | — | 4.363 | 1.202 |
| **Forest weighted mean** | Forest | AGB+BGB carbon | inv-var | 1.143 | **3.386** | 0.834 |
| Sothe et al. SC | Soil | SOC 0–1 m | per-pixel | 35.717 | 42.953 | 27.394 |
| SoilGrids OCS | Soil | SOC 0–1 m | derived | 35.116 | 33.860 | 7.977 |
| Soil equal-weight mean | Soil | SOC 0–1 m | inter-product SD | 10.499 | 37.995 | 14.810 |
| **Soil weighted mean** | Soil | SOC 0–1 m | inv-var | 24.241 | **38.277** | 16.414 |

† AOI-mean statistics from `bestEffort` region reductions vary slightly
(~±0.2 kg/m²) between runs because Earth Engine pyramids to a coarser scale to
respect the pixel cap; the underlying maps are deterministic (seed 42). The
inverse-variance **forest weighted mean is stable at ≈3.39 kg/m²** (AOI) / 3.58
(forested pixels). For an exact, reduction-independent figure, lower the export
scale or remove `bestEffort` in `printStats`.

*(This table is exported verbatim as `Carbon_Summary_Table_v4_6.csv`.)*

### 5.3 Sample allocation

> **Pending re-run with the Forest/Wetland strata.** §4.6 was updated to stratify by
> **Forest + Wetland (GWL_FCS30)** instead of the ESA broad classes. The allocation
> tables below are from the previous ESA-strata run and will change — re-run
> `[▶ RUN ALL]` and paste the two "NEYMAN ALLOCATION" Console blocks (now labelled
> Forest / Wetland) to refresh. Carbon densities, per-product stats, CV, and the
> power analysis (§4–§5.2) are unaffected by the strata definition.

Uncertainty breaks (3 bins, ESA run): forest p33 = 4.33, p67 = 5.20 kg/m²;
soil p33 = 21.76, p67 = 54.26 kg/m².

**Forest — total allocated = 10 (exact).** `[RESULT: refresh with Forest/Wetland run]`

| Stratum | Land cover | Unc bin | N_h (px) | σ_h | n_h |
|---|---|---|---|---|---|
| — | Forest | 0 (low) | — | — | `[RESULT]` |
| — | Forest | 1 (mid) | — | — | `[RESULT]` |
| — | Forest | 2 (high) | — | — | `[RESULT]` |
| — | Wetland | … | — | — | `[RESULT]` |

**Soil — total allocated = 12 (exact).** `[RESULT: refresh with Forest/Wetland run]`

| Stratum | Land cover | Unc bin | N_h (px) | σ_h | n_h |
|---|---|---|---|---|---|
| — | Forest | 0–2 | — | — | `[RESULT]` |
| — | Wetland | 0–2 | — | — | `[RESULT]` |

**Note.** With Forest + Wetland strata, Neyman now routes high-uncertainty wetland
its proportional share even where wetland area is small. Whether wetland draws ≥1
point depends on its area × σ relative to forest; the Console flags any occupied
stratum that receives 0 points.

---

## 6. Outputs (Drive exports)

Run the pipeline, then open the **Tasks** tab and run each export.

**Primary (weighted ensemble):** `Forest_Weighted_Mean_kgm2`,
`Forest_Weighted_Sigma_kgm2`, `Soil_Weighted_Mean_kgm2`, `Soil_Weighted_Sigma_kgm2`,
`SoilGrids_Derived_Sigma_kgm2`, `Total_Ecosystem_Carbon_Weighted_kgm2`.

**Comparison (equal weight):** `Forest_EqualWeight_Mean_kgm2`,
`Forest_EqualWeight_SD_kgm2`, `Soil_EqualWeight_Mean_kgm2`,
`Soil_InterProduct_Disagree_SD_kgm2`, `Total_Ecosystem_Carbon_EqualWeight_kgm2`.

**Uncertainty:** `Sothe_FC_Uncertainty_kgm2`, `Sothe_SC_Uncertainty_kgm2`,
`Forest_RSS_Uncertainty_kgm2`, `Soil_RSS_Uncertainty_kgm2`.

**Source products:** `GEDI_RF_Forest_Carbon_kgm2`, `Sothe_FC_kgm2`,
`SCANFI_Forest_Carbon_kgm2`, `Sothe_SC_kgm2`, `SoilGrids_OCS_0_100cm_kgm2`.

**Data-source maps:** `Forest_Data_Source_Map`, `Soil_Data_Source_Map`.

**Tables / points:** `Carbon_Summary_Table_v4_6.csv`;
`Forest_Sampling_Points_Neyman` and `Soil_Sampling_Points_Neyman` as **KML**
(Google Earth) and **CSV** (lat/lon + carbon, σ, land cover, stratum per point).

---

## 7. How to run (partner instructions)

1. Open `GEE_Script.js` in the [Earth Engine Code Editor](https://code.earthengine.google.com/).
2. In **SECTION A — CONFIGURATION**, set `AOI_ASSET` to your boundary and adjust
   `N_FOREST_SAMPLES` / `N_SOIL_SAMPLES` if needed (these are the **exact** point
   counts produced).
3. Click **`[0] Check Asset Access`**. Every REQUIRED asset must report `✓ PASS`.
   A `✗ FAIL` means that asset must be shared with your Earth Engine account (or its
   ID replaced). Required: the AOI and the two Sothe uncertainty layers.
4. Click **`[▶ RUN ALL]`** and watch the **PROGRESS** panel fill in (Steps 2→9).
   This takes a few minutes.
5. Open the **Tasks** tab (top-right) and click **Run** on each export to save
   rasters, the summary table, and the sampling points to your Google Drive.
6. Load the `*_Sampling_Points_Neyman.kml` files in Google Earth (or the CSVs in any
   GIS) to navigate to each field location.

*(Step `[1]` field-data import and `[4b]` covariate snapshot are optional and are not
needed for the carbon maps or sampling design.)*

---

## 8. Quality assurance & limitations

- **Units harmonised (v4.6).** All forest members are AGB+BGB carbon; previously
  GEDI was biomass and SCANFI was AGB-only carbon, which biased the forest ensemble.
- **Ensemble masking (v4.6).** The weighted mean now honours the "≥2 members"
  rule rather than silently requiring all three present.
- **Exact sample counts (v4.6).** Largest-remainder allocation guarantees the totals
  equal the configured N.
- **BGB approximation.** A single national ratio (0.22) is used for GEDI/SCANFI;
  Sothe used type-specific ratios and also includes dead-plant carbon — a declared,
  modest scope difference.
- **GEDI weighting caveat (material in this AOI).** The GEDI RF σ (3-fold CV RMSE,
  1.358 kg/m²) captures only RF *model* error against GEDI — it excludes GEDI's own
  footprint-level error. Because it is the smallest of the three forest σ values
  (vs SCANFI 2.361 and Sothe 4.721), GEDI RF receives the largest inverse-variance
  weight and **pulls the forest weighted mean (3.39 kg/m²) down toward GEDI's own
  estimate (3.09)** and well below Sothe (5.98). If GEDI's intrinsic uncertainty
  were added to its σ, the weighted mean would shift upward. Field plots will
  resolve which prior is closest. Treat the forest weighted mean as GEDI-leaning.
- **Small forest training set.** Only 110 GEDI L4A footprints fell in the AOI buffer;
  the RF is correspondingly data-limited (CV RMSE ≈ 63% of mean AGBD).
- **Optional asset missing.** `…/combined_profiles` was not found at run time; this
  affects only the Step 4b covariate-extraction CSV, not the carbon maps or sampling.
- **`bestEffort` reductions.** Reported AOI-mean stats can vary ~±0.2 kg/m² between
  runs (Earth Engine pyramids to a coarser scale under the pixel cap); maps and
  exports are deterministic. Lower the scale / drop `bestEffort` for exact figures.
- **Sampling spread is AOI-limited.** With a forest-dominated AOI and a 10/12 budget,
  points concentrate in forest land cover regardless of binning; the design spreads
  them across uncertainty levels rather than land-cover classes.
- **SoilGrids σ is derived**, not a published per-pixel product — treat as indicative.
- **RSS uncertainty (stratification signal)** mixes an absolute SD with a relative
  GEDI term; it is used only to *rank* pixels for stratification, not as the reported
  uncertainty (which is the inverse-variance σ).
- **Power analysis is an SRS bound** using map spatial SD as a proxy for plot
  variance; stratified field sampling should achieve at least this precision.
- **Resolution mixing.** Forest (25 m) and soil (250 m) layers are combined for the
  total; Earth Engine resamples on the fly.

---

## 9. References

- Sothe, C., et al. (2022). *Large Soil Carbon Storage in Terrestrial Ecosystems of
  Canada.* Global Biogeochemical Cycles. https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2021GB007213
- WWF-Canada / McMaster Remote Sensing Lab — Carbon Map. https://wwf.ca/carbonmap/
- Dubayah, R., et al. GEDI L4A Aboveground Biomass Density (ORNL DAAC / LARSE).
- Poggio, L., et al. (2021). *SoilGrids 2.0.* SOIL.
- Guindon, L., et al. (2024). SCANFI: Spatialized CAnadian National Forest Inventory.
- IPCC (2006/2019) Guidelines — biomass carbon fraction and root-to-shoot ratios.
- ESA WorldCover v200 (2021).

---

*Pipeline source: `GEE_Script.js` (v4.6). Results populated from the v4.6 / 3-bin run.
To refresh after any config change, re-run `[▶ RUN ALL]` and update §5 from the
Console (model-performance, total-stocks, and Neyman-allocation blocks).*
