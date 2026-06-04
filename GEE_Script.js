// =================================================================
// Charlie's Place KBA — Forest Carbon Assessment v4.6
// =================================================================
//
//  ┌─────────────────────────────────────────────────────────────┐
//  │  PARTNER SETUP — READ BEFORE RUNNING                          │
//  ├─────────────────────────────────────────────────────────────┤
//  │ 1. ACCESS: this script reads several Earth Engine assets that │
//  │    must be shared with the account running it (or made        │
//  │    public). Click [0] Check Asset Access first — it reports    │
//  │    which assets are readable before you run anything.          │
//  │    Required for the core analysis (Steps 2-9):                 │
//  │      • AOI_ASSET            (your project boundary)            │
//  │      • SOTHE_FC_UNC_ASSET   (Sothe forest-carbon uncertainty)  │
//  │      • SOTHE_SC_UNC_ASSET   (Sothe soil-carbon uncertainty)    │
//  │    Optional (Step 1 / Step 4b field-data exports only):        │
//  │      • WOSIS / CANPEAT / COMBINED assets                       │
//  │ 2. CONFIGURE: set AOI_ASSET to your boundary and adjust        │
//  │    N_FOREST_SAMPLES / N_SOIL_SAMPLES in SECTION A. Whatever    │
//  │    you enter is the EXACT number of points produced.          │
//  │ 3. RUN: click [▶ RUN ALL] to execute Steps 2-9 in order, or    │
//  │    run the numbered buttons one at a time.                     │
//  │ 4. EXPORT: open the Tasks tab (top-right) and click Run on     │
//  │    each export to save rasters / sampling points to Drive.     │
//  └─────────────────────────────────────────────────────────────┘
//
// Changes from v4.5 (v4.6):
//   CARBON UNITS — all forest members harmonised to AGB+BGB carbon
//     (the quantity the Sothe 2022 map represents). New constants
//     BGB_RATIO=0.22, CARBON_FRACTION=0.5, MGHA_TO_KGM2=0.1 →
//     AGB_TO_C_KGM2 = 0.061. GEDI RF and SCANFI (and their σ) now use
//     this factor; previously GEDI was biomass (×0.1, no C fraction)
//     and SCANFI was aboveground-only carbon (×0.05) — incomparable.
//   ENSEMBLE — equal-weight forest ensemble now uses the same 3
//     members as the weighted one (GEDI RF + Sothe FC + SCANFI; was
//     SBFI). Weighted mean now honours the ≥2-members rule via per-
//     member zero-weighting instead of silently requiring all three.
//   dsm_forest now reports # of forest products present (1–3).
//   SAMPLING — Neyman allocation uses largest-remainder rounding so
//     the total ALWAYS equals N_FOREST_SAMPLES / N_SOIL_SAMPLES
//     exactly. Power analysis prints the recommended n vs configured n.
//   PACKAGING — [0] asset-access preflight + [▶ RUN ALL] orchestration
//     so a partner can run end-to-end without ordering steps by hand.
//
// Changes from v4.4:
//
//   Global constants:
//     SCANFI_SIGMA_FC = 1.940 kg C/m²  (Guindon 2024, AGB RMSE
//       38.70 t/ha × 0.1 unit × 0.5 carbon fraction)
//
//   New globals:
//     GEDI_RF_RMSE_MGHA / GEDI_RF_SIGMA_KGPM2 — 3-fold CV output
//     forest_wmean / forest_wmean_sigma — inv-var weighted forest
//     soil_wmean   / soil_wmean_sigma   — inv-var weighted soil
//     sg_sigma     — derived SoilGrids uncertainty image
//     total_ecosystem_c_w — weighted total ecosystem carbon
//
//   Step 5 — 3-fold cross-validation added after final RF training:
//     randomColumn (seed 42) assigns fold 0/1/2 server-side.
//     For each fold: train RF on 2 folds, predict on hold-out,
//     compute fold RMSE. Mean of 3 fold RMSEs → GEDI_RF_RMSE_MGHA.
//     σ in carbon units: GEDI_RF_SIGMA_KGPM2 = GEDI_RF_RMSE_MGHA × 0.1.
//     Step 7 guards on GEDI_RF_SIGMA_KGPM2 being set.
//
//   Step 7 — inverse-variance weighted ensembles replace equal-weight:
//
//     Forest (3 members):
//       GEDI RF  — σ = GEDI_RF_SIGMA_KGPM2 (constant, from 3-fold CV)
//       Sothe FC — σ = sothe_fc_unc (per-pixel kg/m²)
//       SCANFI   — σ = SCANFI_SIGMA_FC (constant 1.940 kg/m²)
//       w_i(p) = 1 / σ_i(p)²
//       forest_wmean = Σ[w_i × x_i] / Σ[w_i]
//       forest_wmean_sigma = 1 / sqrt(Σ[w_i])
//
//     Soil (2 members):
//       Sothe SC — σ = sothe_sc_unc (per-pixel kg/m²)
//       SoilGrids — σ derived (no per-pixel layer available):
//         1. Normalise sothe_sc_unc by its AOI mean → index A (mean=1)
//         2. Normalise inter-product SD (soil_ens_sd) by its AOI mean → index B (mean=1)
//         3. Equal-weight average: combined_index = (A + B) / 2
//         4. Scale by Sothe mean σ → sg_sigma (physical kg/m²)
//         Rationale: SoilGrids is most uncertain where Sothe is uncertain
//         AND where the two products disagree; standardisation ensures
//         neither component dominates by magnitude.
//       soil_wmean = Σ[w_i × x_i] / Σ[w_i]
//       soil_wmean_sigma = 1 / sqrt(Σ[w_i])
//
//     Equal-weight ensemble means retained globally for comparison.
//     forest_nf_mask updated to use forest_wmean.
//     total_ecosystem_c_w = forest_wmean + soil_wmean.
//
//   Step 8 — Neyman allocation uses forest_wmean / soil_wmean.
//             Power analysis uses weighted-mean spatial SD.
//
//   Step 9 — Model performance table: per-product AOI mean, SD,
//             uncertainty metric type and σ value.
//             Total carbon stocks in tonnes and kg added.
//             OOB/CV RMSE interpretation note in console.
// =================================================================


// ─────────────────────────────────────────────────────────────────
// SECTION A — CONFIGURATION
// ─────────────────────────────────────────────────────────────────
var AOI_ASSET      = 'projects/blue-carbon-hub/assets/Charlies_Place_KBA_boundaryFile_2025';
var WOSIS_ASSET    = 'projects/north-star-project-470316/assets/wosis_layers_canada';
var CANPEAT_ASSET  = 'projects/north-star-project-470316/assets/peat_profiles';
var COMBINED_ASSET = 'projects/north-star-project-470316/assets/combined_profiles';

var SOTHE_FC_UNC_ASSET = 'projects/carbon-learning-library/assets/McMasterWWFCanadaforestcarbon250mkgm2uncertaintyversion1';
var SOTHE_SC_UNC_ASSET = 'projects/carbon-learning-library/assets/McMasterWWFCanadasoilcarbon1muncertainty250mkgm2version3';

var EXPORT_CRS     = 'EPSG:32621';
var EXPORT_SCALE   = 25;
var EXPORT_FOLDER  = 'CharliesPlace_Carbon_2025';
var SNAPSHOT_SCALE = 30;

var EMBEDDING_YEAR = null;
var EMBED_SCALE    = 10;

var N_FOREST_SAMPLES = 10;
var N_SOIL_SAMPLES   = 12;
var FOREST_THRESHOLD = 0.5;

var GEDI_START = '2019-04-18';
var GEDI_END   = '2023-12-31';
var S2_START   = '2020-01-01';
var S2_END     = '2023-12-31';
var SAR_START  = '2020-01-01';
var SAR_END    = '2023-12-31';
var TC_START   = '2000-01-01';
var TC_END     = '2022-12-31';

var TOP_N_BANDS       = 6;
var PILOT_N_TREES     = 100;
var FINAL_N_TREES     = 100;
var MIN_GEDI_POINTS   = 50;
var GEDI_SAMPLE_SCALE = 100;
var GEDI_NUM_PIXELS   = 3000;

var N_UNC_BINS          = 4;
var MIN_PTS_PER_STRATUM = 2;   // advisory floor; exact total (above) always wins (v4.6)

// ── Carbon conversion constants (v4.6) ───────────────────────────
// All forest members are harmonised to the SAME quantity as the Sothe et al.
// (2022) forest carbon map: TOTAL biomass carbon = aboveground + belowground.
//   BGB (belowground/root biomass) = BGB_RATIO × AGB
//   Carbon (kg C/m²) = AGB(Mg/ha) × (1 + BGB_RATIO) × CARBON_FRACTION × MGHA_TO_KGM2
// GEDI L4A AGBD and SCANFI biomass are AGB-only, so both receive this factor.
// NOTE: Sothe additionally includes dead-plant carbon and used forest-type-
// specific root ratios; BGB_RATIO here is a single national approximation.
var BGB_RATIO       = 0.22;    // belowground:aboveground biomass ratio
var CARBON_FRACTION = 0.50;    // IPCC dry-biomass carbon fraction
var MGHA_TO_KGM2    = 0.10;    // Mg/ha → kg/m²
var AGB_TO_C_KGM2   = (1 + BGB_RATIO) * CARBON_FRACTION * MGHA_TO_KGM2;  // = 0.061

// SCANFI constant σ = Guindon 2024 AGB RMSE 38.70 t/ha, converted with the
// SAME AGB→carbon factor so it lives in the ensemble's carbon units.
var SCANFI_RMSE_AGB_MGHA = 38.70;                              // Guindon 2024 AGB RMSE (t/ha)
var SCANFI_SIGMA_FC      = SCANFI_RMSE_AGB_MGHA * AGB_TO_C_KGM2;  // ≈ 2.361 kg C/m²

var AGBD_VIS = {
  min: 0, max: 250,
  palette: ['white', 'yellow', 'orange', 'green', 'darkgreen']
};


// ─────────────────────────────────────────────────────────────────
// SECTION B — AOI + CANADA BOUNDARY
// ─────────────────────────────────────────────────────────────────
var aoiFC      = ee.FeatureCollection(AOI_ASSET);
var aoi        = aoiFC.union().geometry();
var aoi_buffer = aoi.buffer(50000);

var canada_boundary = ee.FeatureCollection('USDOS/LSIB_SIMPLE/2017')
  .filter(ee.Filter.eq('country_na', 'Canada'));

Map.centerObject(aoi, 11);
Map.addLayer(canada_boundary, { color: 'b0b0b0', fillColor: '00000000', width: 1 },
  'Canada Boundary', false);
Map.addLayer(aoiFC, { color: '000000', fillColor: '00000000', width: 2 },
  'Conservation Project Boundary', false);


// ─────────────────────────────────────────────────────────────────
// SECTION C — GLOBAL PIPELINE STATE
// ─────────────────────────────────────────────────────────────────
var WOSIS_2023_Raw      = null;
var CanPeatData_Raw     = null;
var sothe_fc            = null;
var sothe_fc_unc        = null;
var sothe_sc            = null;
var sothe_sc_unc        = null;
var scanfi_fc           = null;   // SCANFI FC (kg/m²): AGB × AGB_TO_C_KGM2 (AGB+BGB carbon)
var sg_soc_1m           = null;
var soil_prior_mean     = null;
var scanfi_img          = null;
var gedi_l2a            = null;
var gedi_l4a_col        = null;
var gedi_agbd_se        = null;
var sothe_ch            = null;
var meta_ch             = null;
var scanfi_ch           = null;
var covariates          = null;
var covariates_filled   = null;
var COV_BANDS           = null;
var COV_BANDS_CLIENT    = null;
var HEGL_BANDS          = null;
var sbfi_agb_raster     = null;
var forest_rf_pred      = null;
var forest_ens_mean     = null;   // equal-weight mean (retained for comparison)
var forest_ens_sd       = null;
var soil_ens_mean       = null;   // equal-weight mean (retained for comparison)
var soil_ens_sd         = null;
var forest_uncertainty  = null;   // RSS uncertainty (kept for Neyman reference)
var soil_uncertainty    = null;   // RSS uncertainty (kept for Neyman reference)
var total_ecosystem_c   = null;   // equal-weight total (comparison only)
var forest_nf_mask      = null;   // now based on forest_wmean
var dsm_forest          = null;
var dsm_soil            = null;
var forest_sampling_pts = null;
var soil_sampling_pts   = null;

// 2-stage RF globals
var TOP_BANDS_CLIENT = null;
var TOP_BANDS_EE     = null;
var gedi_training    = null;
var agbd_modeled     = null;

// v4.5: 3-fold CV and inverse-variance weighted ensemble globals
var GEDI_RF_RMSE_MGHA   = null;   // mean 3-fold CV RMSE (Mg/ha AGBD)
var GEDI_RF_SIGMA_KGPM2 = null;   // converted to FC kg/m²: RMSE × AGB_TO_C_KGM2 (AGB+BGB)
var forest_wmean        = null;   // inv-var weighted forest carbon (kg/m²)
var forest_wmean_sigma  = null;   // combined σ from inv-var weighting (kg/m²)
var soil_wmean          = null;   // inv-var weighted soil carbon (kg/m²)
var soil_wmean_sigma    = null;   // combined σ from inv-var weighting (kg/m²)
var sg_sigma            = null;   // derived SoilGrids uncertainty image (kg/m²)
var total_ecosystem_c_w = null;   // forest_wmean + soil_wmean

// snapshot / embedding globals
var snapshot_btn  = null;
var sg_ocs_layers = null;
var embed_img     = null;


// ─────────────────────────────────────────────────────────────────
// SECTION D — UTILITY FUNCTIONS
// ─────────────────────────────────────────────────────────────────
function printStats(img, band, geom, scale, label) {
  img.select(band).reduceRegion({
    reducer: ee.Reducer.mean()
      .combine(ee.Reducer.min(),    null, true)
      .combine(ee.Reducer.max(),    null, true)
      .combine(ee.Reducer.stdDev(), null, true),
    geometry: geom  || aoi,
    scale:    scale || EXPORT_SCALE,
    crs:      EXPORT_CRS,
    maxPixels: 1e11, bestEffort: true
  }).evaluate(function(r) {
    var f = function(k) { return r[k] !== undefined ? r[k].toFixed(3) : 'null'; };
    print(label + ' — mean: ' + f(band + '_mean') +
          ' | min: '  + f(band + '_min') +
          ' | max: '  + f(band + '_max') +
          ' | sd: '   + f(band + '_stdDev'));
  });
}

function addLayerWithStats(img, band, vis, label, geom, scale) {
  Map.addLayer(img.select(band), vis, label, false);
  printStats(img, band, geom, scale, label);
}

function getStats(img, band, sc) {
  var s = img.select(band).reduceRegion({
    reducer: ee.Reducer.mean()
      .combine(ee.Reducer.stdDev(), null, true)
      .combine(ee.Reducer.min(),    null, true)
      .combine(ee.Reducer.max(),    null, true),
    geometry: aoi,
    scale:    sc || EXPORT_SCALE,
    crs:      EXPORT_CRS,
    maxPixels: 1e11, bestEffort: true
  });
  return {
    mean: s.get(band + '_mean'),
    sd:   s.get(band + '_stdDev'),
    min:  s.get(band + '_min'),
    max:  s.get(band + '_max')
  };
}

function cvStat(sd, mean) {
  return ee.Number(sd).divide(ee.Number(mean).abs().add(0.001)).multiply(100);
}


// ─────────────────────────────────────────────────────────────────
// SECTION E — SAR UTILITY
// ─────────────────────────────────────────────────────────────────
function buildSARComposites() {
  var s1 = ee.ImageCollection('COPERNICUS/S1_GRD')
    .filterBounds(aoi_buffer).filterDate(SAR_START, SAR_END)
    .filter(ee.Filter.eq('instrumentMode', 'IW'))
    .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'))
    .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH'))
    .select(['VV', 'VH'])
    .map(function(img) { return img.updateMask(img.select('VV').gt(-30)); });

  var asc  = s1.filter(ee.Filter.eq('orbitProperties_pass', 'ASCENDING'));
  var desc = s1.filter(ee.Filter.eq('orbitProperties_pass', 'DESCENDING'));
  var emptyOrbit = ee.Image.constant([-999,-999]).rename(['VV','VH']).updateMask(ee.Image(0));

  var seasonal = function(col, m0, m1) {
    var sub = col.filter(ee.Filter.calendarRange(m0, m1, 'month'));
    return ee.Image(ee.Algorithms.If(sub.size().gt(0), sub.median(), emptyOrbit));
  };

  var asc_spring  = seasonal(asc,  3, 4);
  var asc_summer  = seasonal(asc,  7, 8);
  var desc_spring = seasonal(desc, 3, 4);
  var desc_summer = seasonal(desc, 7, 8);

  return {
    asc_spring:  asc_spring,  asc_summer:  asc_summer,
    desc_spring: desc_spring, desc_summer: desc_summer,
    asc_change:  asc_summer.select('VV').subtract(asc_spring.select('VV')).rename('asc_vv_change'),
    desc_change: desc_summer.select('VV').subtract(desc_spring.select('VV')).rename('desc_vv_change')
  };
}


// ─────────────────────────────────────────────────────────────────
// SECTION F — UI
// ─────────────────────────────────────────────────────────────────
var statusLabel    = null;
var progressLabels = [];
var progressSteps  = [
  '[ ] Step 1: Raw data imported',
  '[ ] Step 2: Raster priors loaded',
  '[ ] Step 3: GEDI + CHMs + SCANFI FC',
  '[ ] Step 4: Covariates + variable selection',
  '[ ] Step 5: Forest model + 3-fold CV',
  '[ ] Step 6: Soil sources verified',
  '[ ] Step 7: Inv-var weighted ensembles',
  '[ ] Step 8: Sampling complete',
  '[ ] Step 9: Exports queued'
];
var mainPanel = null;

function markDone(idx) {
  if (!progressLabels[idx]) return;
  progressLabels[idx].setValue(progressSteps[idx].replace('[ ]', '[✓]'));
  progressLabels[idx].style().set('color', '#2e7d32');
}

function setStatus(msg) {
  if (statusLabel) statusLabel.setValue(msg);
  print('STATUS: ' + msg);
}

function makeButton(label, onClick) {
  return ui.Button({ label: label, onClick: onClick,
    style: { width: '260px', margin: '2px 0', fontSize: '11px' } });
}

function makeSectionLabel(text) {
  return ui.Label(text, {
    fontWeight: 'bold', fontSize: '11px', color: '#222222', margin: '8px 0 3px 0'
  });
}

var SNAPSHOT_BTN_INDEX = 7;

function insertSnapshotButton() {
  if (snapshot_btn) return;
  snapshot_btn = ui.Button({
    label: '[4b] Export Covariate Snapshot (30 m GeoTIFF)',
    onClick: step4b_exportCovariateSnapshot,
    style: { width: '260px', margin: '2px 0', fontSize: '11px',
             color: '#ffffff', backgroundColor: '#1565c0' }
  });
  if (mainPanel) { mainPanel.insert(SNAPSHOT_BTN_INDEX, snapshot_btn); }
}

function initUI() {
  mainPanel = ui.Panel({
    style: { width: '290px', padding: '14px', position: 'top-left', backgroundColor: '#ffffff' }
  });
  mainPanel.add(ui.Label("Charlie's Place KBA", {
    fontWeight: 'bold', fontSize: '16px', color: '#111111', margin: '0 0 2px 0'
  }));
  mainPanel.add(ui.Label('Forest Carbon Assessment v4.6', {
    fontSize: '12px', color: '#555555', margin: '0 0 12px 0'
  }));

  mainPanel.add(makeSectionLabel('SETUP / ONE-CLICK'));
  mainPanel.add(makeButton('[0] Check Asset Access', checkAssetAccess));
  mainPanel.add(ui.Button({
    label: '[▶ RUN ALL]  Steps 2 → 9',
    onClick: runAll,
    style: { width: '260px', margin: '2px 0 8px 0', fontSize: '11px',
             color: '#ffffff', backgroundColor: '#2e7d32', fontWeight: 'bold' }
  }));

  mainPanel.add(makeSectionLabel('PIPELINE (manual / step-by-step)'));
  mainPanel.add(makeButton('[1] Import Raw Field Data (optional)',  step1_importRawData));
  mainPanel.add(makeButton('[2] Import Raster Priors',              step2_importRasterPriors));
  mainPanel.add(makeButton('[3] GEDI + CHMs + SCANFI FC',           step3_importGEDI));
  mainPanel.add(makeButton('[4] Build Covariates + Var Selection',   step4_buildCovariates));
  mainPanel.add(makeButton('[5] Train Final RF + 3-fold CV',        step5_trainFinalGEDI));
  mainPanel.add(makeButton('[6] Verify Soil Sources',               step6_buildSoilModel));
  mainPanel.add(makeButton('[7] Build Inv-Var Weighted Ensembles',  step7_buildEnsemble));
  mainPanel.add(makeButton('[8] Generate Sampling Points',          step8_generateSampling));
  mainPanel.add(makeButton('[9] Reports and Exports',               step9_reportsAndExports));

  mainPanel.add(makeSectionLabel('PROGRESS'));
  var progressPanel = ui.Panel({ style: { margin: '0 0 8px 0' } });
  progressLabels = progressSteps.map(function(txt) {
    var lbl = ui.Label(txt, { fontSize: '10px', color: '#999999', margin: '1px 0' });
    progressPanel.add(lbl);
    return lbl;
  });
  mainPanel.add(progressPanel);
  mainPanel.add(makeSectionLabel('STATUS'));
  statusLabel = ui.Label('Ready - click [1] to begin.', {
    fontSize: '11px', color: '#333333', margin: '0 0 10px 0'
  });
  mainPanel.add(statusLabel);
  mainPanel.add(ui.Label(
    'Scale: ' + EXPORT_SCALE + 'm | CRS: UTM Zone 21N\n' +
    'Forest: ' + N_FOREST_SAMPLES + ' pts | Soil: ' + N_SOIL_SAMPLES + ' pts (exact)\n' +
    'Carbon: AGB+BGB | AGB→C ×' + AGB_TO_C_KGM2.toFixed(3) + '\n' +
    'Ensemble: inv-var weighting (1/σ²) per pool\n' +
    'Sampling: Neyman (exact n) | LC: ESA WorldCover | unc bins: ' + N_UNC_BINS,
    { fontSize: '9px', color: '#aaaaaa', margin: '4px 0 0 0', whiteSpace: 'pre' }
  ));
  Map.add(mainPanel);
}


// ─────────────────────────────────────────────────────────────────
// STEP 4b — EXPORT COVARIATE SNAPSHOT
// ─────────────────────────────────────────────────────────────────
function step4b_exportCovariateSnapshot() {
  if (!covariates)   { setStatus('Run Step 4 first.'); return; }
  if (!sg_ocs_layers){ setStatus('Run Step 2 first.'); return; }
  if (!sothe_sc)     { setStatus('Run Step 2 first.'); return; }
  if (!sothe_fc)     { setStatus('Run Step 2 first.'); return; }
  if (!sothe_ch)     { setStatus('Run Step 3 first.'); return; }
  if (!gedi_l4a_col) { setStatus('Run Step 3 first.'); return; }
  if (!embed_img)    { setStatus('Run Step 2 first.'); return; }

  setStatus('Step 4b: Assembling covariate snapshot stack...');

  var gedi_l4a_median = gedi_l4a_col
    .map(function(img) {
      var q = img.select('l4_quality_flag').eq(1).and(img.select('degrade_flag').eq(0));
      return img.select(['agbd', 'agbd_se']).updateMask(q);
    }).median().rename(['gedi_agbd', 'agbd_se']);

  var cov_stack_global = covariates
    .addBands(sg_ocs_layers)
    .addBands(sothe_sc.rename('sothe_sc'))
    .addBands(sothe_fc.rename('sothe_fc'))
    .addBands(sothe_ch.select([0]).rename('sothe_ch'))
    .addBands(gedi_l4a_median);
  if (scanfi_fc) cov_stack_global = cov_stack_global.addBands(scanfi_fc.rename('scanfi_fc'));

  var snapshot = cov_stack_global.clip(aoi);

  snapshot.bandNames().evaluate(function(names) {
    print('COVARIATE SNAPSHOT BAND INVENTORY (' + names.length + ' bands)');
    names.forEach(function(b, i) { print('  ' + (i+1) + '. ' + b); });
  });

  Export.image.toDrive({ image: snapshot.toFloat(), description: 'CharliesPlace_Covariate_Snapshot_30m',
    folder: EXPORT_FOLDER, region: aoi, scale: SNAPSHOT_SCALE, crs: EXPORT_CRS, maxPixels: 1e13 });
  Export.image.toDrive({ image: embed_img.clip(aoi).toFloat(), description: 'CharliesPlace_GoogleEmbedding_V1_10m',
    folder: EXPORT_FOLDER, region: aoi, scale: EMBED_SCALE, crs: EXPORT_CRS, maxPixels: 1e13 });

  var combined     = ee.FeatureCollection(COMBINED_ASSET);
  var wosis_na     = combined.filter(ee.Filter.eq('dataset', 'WOSIS 2023'))
                             .filter(ee.Filter.inList('country_name', ee.List(['Canada','United States of America'])));
  var canpeat_pts  = combined.filter(ee.Filter.eq('dataset', 'Peat Database'));
  var janousek_pts = combined.filter(ee.Filter.eq('dataset', 'Janousek'));
  var core_pts     = wosis_na.merge(canpeat_pts).merge(janousek_pts);

  Export.table.toDrive({ collection: cov_stack_global.sampleRegions({
    collection: core_pts, scale: SNAPSHOT_SCALE, tileScale: 4, geometries: true }),
    description: 'CorePoints_Covariates_CSV', folder: EXPORT_FOLDER, fileFormat: 'CSV' });
  Export.table.toDrive({ collection: embed_img.sampleRegions({
    collection: core_pts, scale: EMBED_SCALE, tileScale: 4, geometries: true }),
    description: 'CorePoints_GoogleEmbedding_V1_CSV', folder: EXPORT_FOLDER, fileFormat: 'CSV' });

  setStatus('Step 4b: 4 tasks queued. Open Tasks panel to run.');
}


// ─────────────────────────────────────────────────────────────────
// STEP 1 — IMPORT RAW FIELD DATA
// ─────────────────────────────────────────────────────────────────
function step1_importRawData() {
  setStatus('Step 1: Loading WOSIS and CanPeat field data...');
  WOSIS_2023_Raw  = ee.FeatureCollection(WOSIS_ASSET);
  CanPeatData_Raw = ee.FeatureCollection(CANPEAT_ASSET);
  WOSIS_2023_Raw.size().evaluate(function(n)  { print('WOSIS_2023_Raw - count:', n); });
  CanPeatData_Raw.size().evaluate(function(n) { print('CanPeatData_Raw - count:', n); });
  Map.addLayer(WOSIS_2023_Raw,  { color: '1565c0' }, 'WOSIS_2023_Raw',  false);
  Map.addLayer(CanPeatData_Raw, { color: 'e65100' }, 'CanPeatData_Raw', false);
  markDone(0);
  setStatus('Step 1 complete.');
}


// ─────────────────────────────────────────────────────────────────
// STEP 2 — RASTER PRIORS
// ─────────────────────────────────────────────────────────────────
function step2_importRasterPriors(onDone) {
  setStatus('Step 2: Loading raster priors...');

  sothe_fc = ee.ImageCollection('projects/sat-io/open-datasets/carbon_stocks_ca/fc')
    .first().clip(aoi).rename('sothe_fc');
  sothe_sc = ee.ImageCollection('projects/sat-io/open-datasets/carbon_stocks_ca/sc')
    .first().clip(aoi).rename('sothe_sc');
  sothe_fc.bandNames().evaluate(function(n) { print('Sothe FC bands:', n); });
  sothe_sc.bandNames().evaluate(function(n) { print('Sothe SC bands:', n); });

  sothe_fc_unc = ee.Image(SOTHE_FC_UNC_ASSET).select([0]).clip(aoi).rename('sothe_fc_unc');
  sothe_sc_unc = ee.Image(SOTHE_SC_UNC_ASSET).select([0]).clip(aoi).rename('sothe_sc_unc');
  sothe_fc_unc.bandNames().evaluate(function(n) { print('Sothe FC unc band:', n); });
  sothe_sc_unc.bandNames().evaluate(function(n) { print('Sothe SC unc band:', n); });

  var sg_soc  = ee.Image('projects/soilgrids-isric/soc_mean').clip(aoi_buffer);
  var sg_bdod = ee.Image('projects/soilgrids-isric/bdod_mean').clip(aoi_buffer);
  sg_soc.bandNames().evaluate(function(n)  { print('SoilGrids SOC bands:',  n); });
  sg_bdod.bandNames().evaluate(function(n) { print('SoilGrids BDOD bands:', n); });

  var ocsLayer = function(socBand, bdBand, thickness) {
    return sg_soc.select(socBand).divide(10)
      .multiply(sg_bdod.select(bdBand).divide(100))
      .multiply(thickness).divide(100);
  };
  var ocs_0_5    = ocsLayer('soc_0-5cm_mean',    'bdod_0-5cm_mean',     5).rename('sg_ocs_0_5cm');
  var ocs_5_15   = ocsLayer('soc_5-15cm_mean',   'bdod_5-15cm_mean',   10).rename('sg_ocs_5_15cm');
  var ocs_15_30  = ocsLayer('soc_15-30cm_mean',  'bdod_15-30cm_mean',  15).rename('sg_ocs_15_30cm');
  var ocs_30_60  = ocsLayer('soc_30-60cm_mean',  'bdod_30-60cm_mean',  30).rename('sg_ocs_30_60cm');
  var ocs_60_100 = ocsLayer('soc_60-100cm_mean', 'bdod_60-100cm_mean', 40).rename('sg_ocs_60_100cm');

  sg_ocs_layers = ocs_0_5.addBands(ocs_5_15).addBands(ocs_15_30).addBands(ocs_30_60).addBands(ocs_60_100).clip(aoi);
  sg_soc_1m = ocs_0_5.add(ocs_5_15).add(ocs_15_30).add(ocs_30_60).add(ocs_60_100).rename('sg_soc_1m').clip(aoi);
  soil_prior_mean = sg_soc_1m.rename('sg').addBands(sothe_sc.rename('sothe'))
    .reduce(ee.Reducer.mean()).rename('soil_prior_mean');

  scanfi_img = ee.Image('projects/gcpm041u-lemur/assets/scanfi_v12/SCANFI_v1_2').clip(aoi);
  scanfi_img.bandNames().evaluate(function(n) {
    print('SCANFI v1.2 bands:', n);
    Map.addLayer(scanfi_img.select(n[0]),
      { min: 0, max: 200, palette: ['#f7f7f7', '#74c476', '#00441b'] },
      'SCANFI Biomass (' + n[0] + ') (Mg/ha)', false);
  });

  var soilVis   = { min: 0, max: 30,  palette: ['#fff7bc', '#fe9929', '#993404'] };
  var forestVis = { min: 0, max: 20,  palette: ['#f7fcf5', '#74c476', '#00441b'] };
  var diffVis   = { min: -10, max: 10, palette: ['#d73027', '#ffffbf', '#1a9850'] };
  var uncVis    = { min: 0,  max: 10,  palette: ['#2166ac', '#f7f7f7', '#d6604d'] };

  addLayerWithStats(sg_soc_1m,      'sg_soc_1m',         soilVis,   'SoilGrids SOC 0-100 cm (kg/m2)',          aoi, 250);
  addLayerWithStats(sothe_sc,       'sothe_sc',           soilVis,   'Sothe et al. Soil Carbon (kg/m2)',         aoi, 250);
  addLayerWithStats(sothe_fc,       'sothe_fc',           forestVis, 'Sothe et al. Forest Carbon (kg/m2)',      aoi, 250);
  addLayerWithStats(sothe_fc_unc,   'sothe_fc_unc',       uncVis,    'Sothe FC Uncertainty (kg/m2)',            aoi, 250);
  addLayerWithStats(sothe_sc_unc,   'sothe_sc_unc',       uncVis,    'Sothe SC Uncertainty (kg/m2)',            aoi, 250);
  addLayerWithStats(soil_prior_mean,'soil_prior_mean',    soilVis,   'Soil Prior Mean - equal weight (kg/m2)',  aoi, 250);

  var diff_sc = sothe_sc.subtract(sg_soc_1m).rename('sothe_minus_sg_soc');
  addLayerWithStats(diff_sc, 'sothe_minus_sg_soc', diffVis, 'Soil Diff: Sothe - SoilGrids (kg/m2)', aoi, 250);

  var embedCol = ee.ImageCollection('GOOGLE/SATELLITE_EMBEDDING/V1/ANNUAL').filterBounds(aoi_buffer);
  if (EMBEDDING_YEAR !== null) embedCol = embedCol.filter(ee.Filter.calendarRange(EMBEDDING_YEAR, EMBEDDING_YEAR, 'year'));
  embed_img = embedCol.median().clip(aoi_buffer);
  embed_img.bandNames().size().evaluate(function(n) {
    print('Google Satellite Embedding V1 - bands: ' + n + ' | scale: ' + EMBED_SCALE + ' m');
  });

  markDone(1);
  setStatus('Step 2 complete - raster priors, Sothe uncertainty, and Google Embedding loaded.');
  if (typeof onDone === 'function') onDone();
}


// ─────────────────────────────────────────────────────────────────
// STEP 3 — GEDI + CANOPY HEIGHT MODELS + SCANFI FOREST CARBON
//
// SCANFI AGB (band[0], Mg/ha) → forest carbon (kg/m², AGB+BGB):
//   AGB × (1+BGB_RATIO) × CARBON_FRACTION × MGHA_TO_KGM2 = AGB × AGB_TO_C_KGM2 (= 0.061)
//   (v4.6: was AGB × 0.05 — aboveground-only carbon, inconsistent with Sothe)
// ─────────────────────────────────────────────────────────────────
function step3_importGEDI(onDone) {
  setStatus('Step 3: Loading GEDI L2A and L4A...');

  var l2a_col = ee.ImageCollection('LARSE/GEDI/GEDI02_A_002_MONTHLY')
    .filterDate(GEDI_START, GEDI_END).filterBounds(aoi_buffer)
    .map(function(img) {
      var q = img.select('quality_flag').eq(1).and(img.select('degrade_flag').eq(0));
      return img.select(['rh98']).updateMask(q);
    });

  gedi_l4a_col = ee.ImageCollection('LARSE/GEDI/GEDI04_A_002_MONTHLY')
    .filterDate(GEDI_START, GEDI_END).filterBounds(aoi_buffer)
    .select(['agbd', 'agbd_se', 'l4_quality_flag', 'degrade_flag']);

  l2a_col.size().evaluate(function(n)      { print('GEDI L2A images in buffer:', n); });
  gedi_l4a_col.size().evaluate(function(n) { print('GEDI L4A images in buffer:', n); });

  gedi_l2a     = l2a_col.median().clip(aoi);
  // GEDI AGBD SE → carbon σ (AGB+BGB, kg C/m²), same factor as the mean
  gedi_agbd_se = gedi_l4a_col.select('agbd_se').mean().multiply(AGB_TO_C_KGM2).rename('agbd_se_kgm2').clip(aoi);

  var shot_density = l2a_col.map(function(img) {
    return img.select('rh98').mask().rename('shot');
  }).sum().rename('gedi_shot_density').clip(aoi);

  sothe_ch = ee.ImageCollection('projects/sat-io/open-datasets/carbon_stocks_ca/ch')
    .filterBounds(aoi_buffer).first().clip(aoi);
  var sothe_ch_b1 = sothe_ch.select([0]).rename('sothe_ch');

  meta_ch = ee.ImageCollection('projects/sat-io/open-datasets/facebook/meta-canopy-height')
    .filterBounds(aoi_buffer).mosaic()
    .reproject({ crs: EXPORT_CRS, scale: EXPORT_SCALE }).clip(aoi).rename('meta_ch');

  scanfi_ch = scanfi_img.select('height').rename('scanfi_ch');

  printStats(gedi_l2a,   'rh98',      aoi, 25, 'GEDI L2A rh98 Canopy Height (m)');
  printStats(sothe_ch_b1,'sothe_ch',  aoi, 25, 'Sothe Canopy Height Model (m)');
  printStats(meta_ch,    'meta_ch',   aoi, 25, 'Meta Canopy Height Model (m)');
  printStats(scanfi_ch,  'scanfi_ch', aoi, 25, 'SCANFI Canopy Height (m)');

  var gedi_rh98 = gedi_l2a.select('rh98');
  var printMeanDiff = function(chm, chmBand, label) {
    chm.select(chmBand).subtract(gedi_rh98).rename('diff')
      .reduceRegion({ reducer: ee.Reducer.mean(), geometry: aoi,
        scale: EXPORT_SCALE, crs: EXPORT_CRS, maxPixels: 1e11, bestEffort: true })
      .evaluate(function(r) {
        print(label + ' vs GEDI rh98 - mean diff (m): ' +
          (r['diff'] !== undefined ? r['diff'].toFixed(3) : 'null'));
      });
  };
  printMeanDiff(sothe_ch_b1, 'sothe_ch',  'Sothe CHM');
  printMeanDiff(meta_ch,     'meta_ch',   'Meta CHM');
  printMeanDiff(scanfi_ch,   'scanfi_ch', 'SCANFI CHM');

  // SCANFI AGB → Forest Carbon (AGB+BGB, kg C/m²) via AGB_TO_C_KGM2
  scanfi_fc = scanfi_img.select([0]).multiply(AGB_TO_C_KGM2).rename('scanfi_fc').clip(aoi);
  addLayerWithStats(scanfi_fc, 'scanfi_fc',
    { min: 0, max: 20, palette: ['#f7fcf5', '#74c476', '#00441b'] },
    'SCANFI Forest Carbon (kg/m2) [AGB+BGB; σ=' + SCANFI_SIGMA_FC.toFixed(3) + ' constant]', aoi, 25);

  var fc_diff = sothe_fc.subtract(scanfi_fc).rename('sothe_minus_scanfi_fc');
  addLayerWithStats(fc_diff, 'sothe_minus_scanfi_fc',
    { min: -5, max: 5, palette: ['#d73027', '#ffffbf', '#1a9850'] },
    'FC Diff: Sothe - SCANFI (kg/m2)', aoi, 25);

  var htVis = { min: 0, max: 30, palette: ['#f7fcb9', '#addd8e', '#31a354', '#006837'] };
  Map.addLayer(gedi_l2a.select('rh98'), htVis, 'GEDI L2A rh98 (m)', false);
  Map.addLayer(sothe_ch_b1, htVis, 'Sothe CHM (m)',  false);
  Map.addLayer(meta_ch,     htVis, 'Meta CHM (m)',   false);
  Map.addLayer(scanfi_ch,   htVis, 'SCANFI CHM (m)', false);
  Map.addLayer(shot_density, { min: 0, max: 50, palette: ['#ffffcc', '#c7e9b4', '#0c2c84'] },
    'GEDI Shot Density', false);

  markDone(2);
  setStatus('Step 3 complete - GEDI, CHMs, and SCANFI FC loaded.');
  if (typeof onDone === 'function') onDone();
}


// ─────────────────────────────────────────────────────────────────
// STEP 4 — BUILD COVARIATE STACK + PILOT RF
// ─────────────────────────────────────────────────────────────────
function step4_buildCovariates(onDone) {
  if (!gedi_l2a) { setStatus('Run Step 3 first.'); return; }
  setStatus('Step 4: Building covariate stack...');

  var dem    = ee.Image('NASA/NASADEM_HGT/001').select('elevation');
  var slope  = ee.Terrain.slope(dem);
  var aspect = ee.Terrain.aspect(dem);
  var twi    = slope.multiply(Math.PI / 180).tan().max(ee.Image(0.001)).pow(-1).log().rename('twi');
  var tpi    = dem.subtract(dem.reduceNeighborhood({
    reducer: ee.Reducer.mean(), kernel: ee.Kernel.circle({ radius: 500, units: 'meters' })
  })).rename('tpi');

  var s2 = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
    .filterBounds(aoi_buffer).filterDate(S2_START, S2_END)
    .filter(ee.Filter.calendarRange(6, 9, 'month'))
    .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 20))
    .map(function(img) {
      var qa   = img.select('QA60');
      var mask = qa.bitwiseAnd(1 << 10).eq(0).and(qa.bitwiseAnd(1 << 11).eq(0));
      return img.updateMask(mask).divide(10000).copyProperties(img, ['system:time_start']);
    }).median();

  var ndvi = s2.normalizedDifference(['B8', 'B4']).rename('ndvi');
  var ndwi = s2.normalizedDifference(['B3', 'B8']).rename('ndwi');
  var nbr  = s2.normalizedDifference(['B8', 'B12']).rename('nbr');
  var evi  = s2.expression('2.5*((NIR-RED)/(NIR+6*RED-7.5*BLUE+1))',
    { NIR: s2.select('B8'), RED: s2.select('B4'), BLUE: s2.select('B2') }).rename('evi');

  var climate = ee.ImageCollection('IDAHO_EPSCOR/TERRACLIMATE')
    .filterBounds(aoi_buffer).filterDate(TC_START, TC_END)
    .select(['tmmn', 'tmmx', 'pr', 'soil']).mean();

  setStatus('Step 4: Building SAR composites...');
  var sar = buildSARComposites();

  var grid_fe = ee.FeatureCollection(
    'projects/sat-io/open-datasets/CA_FOREST/CA_SBFI/GRID_forested_ecosystems');
  sbfi_agb_raster = grid_fe.filter(ee.Filter.notNull(['STRUCTURE_AGB_AVG']))
    .reduceToImage({ properties: ['STRUCTURE_AGB_AVG'], reducer: ee.Reducer.first() })
    .rename('sbfi_agb_avg')
    .addBands(grid_fe.filter(ee.Filter.notNull(['STRUCTURE_AGB_SD']))
      .reduceToImage({ properties: ['STRUCTURE_AGB_SD'], reducer: ee.Reducer.first() })
      .rename('sbfi_agb_sd'))
    .reproject({ crs: EXPORT_CRS, scale: EXPORT_SCALE }).clip(aoi);
  sbfi_agb_raster.bandNames().evaluate(function(n) { print('SBFI rasterised bands:', n); });
  Map.addLayer(grid_fe, { color: '888888' }, 'SBFI Grid Forested Ecosystems', false);

  var lead_species = ee.Image('projects/sat-io/open-datasets/CA_FOREST/LEAD_TREE_SPECIES')
    .rename('lead_tree_species').clip(aoi_buffer);

  var candidate = dem.rename('elevation')
    .addBands(slope.rename('slope')).addBands(aspect.rename('aspect'))
    .addBands(twi).addBands(tpi)
    .addBands(ndvi).addBands(ndwi).addBands(nbr).addBands(evi)
    .addBands(s2.select(['B2', 'B3', 'B4', 'B8', 'B11', 'B12']))
    .addBands(climate.select('tmmn').rename('tmmn'))
    .addBands(climate.select('tmmx').rename('tmmx'))
    .addBands(climate.select('pr').rename('map'))
    .addBands(climate.select('soil').rename('soil_moisture'))
    .addBands(sar.asc_spring.select('VV').rename('sar_asc_spr_vv'))
    .addBands(sar.asc_spring.select('VH').rename('sar_asc_spr_vh'))
    .addBands(sar.asc_summer.select('VV').rename('sar_asc_sum_vv'))
    .addBands(sar.asc_summer.select('VH').rename('sar_asc_sum_vh'))
    .addBands(sar.asc_change.rename('sar_asc_vv_change'))
    .addBands(sbfi_agb_raster.select('sbfi_agb_avg'))
    .addBands(sbfi_agb_raster.select('sbfi_agb_sd'))
    .addBands(lead_species)
    .reproject({ crs: EXPORT_CRS, scale: EXPORT_SCALE });

  setStatus('Step 4: Checking band coverage over AOI...');

  candidate.reduceRegion({ reducer: ee.Reducer.count(), geometry: aoi,
    scale: EXPORT_SCALE, crs: EXPORT_CRS, maxPixels: 1e11, bestEffort: true
  }).evaluate(function(counts) {
    var allBands = candidate.bandNames().getInfo();
    var goodBands = [], deadBands = [];
    allBands.forEach(function(b) {
      if (counts[b] && counts[b] > 0) goodBands.push(b); else deadBands.push(b);
    });
    print('COVARIATE COVERAGE: ' + goodBands.length + '/' + allBands.length + ' bands with data');
    if (deadBands.length > 0) print('DROPPED bands:', deadBands);

    var coreBands   = ['elevation', 'ndvi', 'B8'];
    var missingCore = coreBands.filter(function(b) { return goodBands.indexOf(b) === -1; });
    if (missingCore.length > 0) { setStatus('Core bands missing: ' + missingCore.join(', ')); return; }

    covariates       = candidate.select(goodBands);
    covariates_filled = covariates.unmask(-9999);
    COV_BANDS        = ee.List(goodBands);
    COV_BANDS_CLIENT = goodBands;
    HEGL_BANDS       = COV_BANDS.cat(ee.List(['sothe_sc_covariate', 'depth_cm']));

    print('Final covariate stack (' + goodBands.length + ' bands):', goodBands);
    covariates.projection().evaluate(function(p) { print('Covariate projection:', p); });

    Map.addLayer(ndvi.clip(aoi), { min: -0.1, max: 0.9, palette: ['#d73027','#ffffbf','#1a9850'] },
      'NDVI (S2 Summer Median)', false);
    Map.addLayer(covariates.select('twi').clip(aoi), { min: 0, max: 10, palette: ['#ffffcc','#41b6c4','#0c2c84'] },
      'TWI', false);
    Map.addLayer(covariates.select('tpi').clip(aoi), { min: -30, max: 30, palette: ['#4575b4','#ffffbf','#d73027'] },
      'TPI', false);

    _step4_pilotRF(goodBands, onDone);
  });
}

function _step4_pilotRF(goodBands, onDone) {
  var gedi_raw = ee.ImageCollection('LARSE/GEDI/GEDI04_A_002_MONTHLY')
    .filterDate(GEDI_START, GEDI_END).filterBounds(aoi_buffer)
    .map(function(img) {
      var q = img.select('l4_quality_flag').eq(1).and(img.select('degrade_flag').eq(0));
      return img.select(['agbd', 'agbd_se']).updateMask(q);
    });

  var sampled = gedi_raw.median().select('agbd').addBands(covariates)
    .sample({ region: aoi_buffer, scale: GEDI_SAMPLE_SCALE, numPixels: GEDI_NUM_PIXELS,
              seed: 42, tileScale: 4, geometries: true, dropNulls: true })
    .filter(ee.Filter.gt('agbd', 0)).filter(ee.Filter.lt('agbd', 600));

  sampled.size().evaluate(function(nPoints) {
    print('GEDI L4A training points (0 < agbd < 600):', nPoints);
    if (nPoints === 0) { setStatus('No GEDI L4A points found.'); return; }

    var trainingSamples;
    if (nPoints < MIN_GEDI_POINTS) {
      print('Low density (' + nPoints + ' pts) - buffering shots by 500 m.');
      var buffered = sampled.map(function(f) { return f.buffer(500); });
      trainingSamples = covariates.select(goodBands).sampleRegions({
        collection: buffered, properties: ['agbd'], scale: EXPORT_SCALE, tileScale: 4
      }).filter(ee.Filter.notNull(goodBands)).filter(ee.Filter.gt('agbd', 0));
    } else {
      trainingSamples = sampled.filter(ee.Filter.notNull(goodBands));
    }
    print('Training samples (complete cases):', trainingSamples.size());
    setStatus('Step 4: Training pilot RF (' + PILOT_N_TREES + ' trees, ' + goodBands.length + ' bands)...');

    var goodBandsEE = ee.List(goodBands);
    var pilotRF = ee.Classifier.smileRandomForest({
      numberOfTrees: PILOT_N_TREES, minLeafPopulation: 5, bagFraction: 0.632, seed: 42
    }).setOutputMode('REGRESSION')
      .train({ features: trainingSamples, classProperty: 'agbd', inputProperties: goodBandsEE });

    var importanceDict = ee.Dictionary(pilotRF.explain().get('importance'));
    ee.FeatureCollection(importanceDict.keys().map(function(key) {
      return ee.Feature(null, { band: key, importance: importanceDict.getNumber(key) });
    })).sort('importance', false).limit(TOP_N_BANDS).aggregate_array('band').evaluate(function(topBands) {
      if (!topBands || topBands.length === 0) { setStatus('Variable selection returned no bands.'); return; }

      importanceDict.evaluate(function(impObj) {
        var chartData = [];
        for (var b in impObj) chartData.push({ band: b, importance: impObj[b] });
        chartData.sort(function(a, b) { return b.importance - a.importance; });
        print(ui.Chart.feature.byFeature(
          ee.FeatureCollection(chartData.map(function(d) { return ee.Feature(null, d); })),
          'band', 'importance').setChartType('BarChart').setOptions({
            title: 'Pilot RF - Variable Importance (all ' + goodBands.length + ' bands)',
            hAxis: { title: 'Importance (Gini)' }, vAxis: { title: 'Band' },
            legend: 'none', colors: ['#2E7D32'], bar: { groupWidth: '80%' }
        }));
      });

      print('Variable selection - top ' + topBands.length + ' bands:');
      topBands.forEach(function(b) { print('  ' + b); });

      TOP_BANDS_CLIENT = topBands;
      TOP_BANDS_EE     = ee.List(topBands);
      gedi_training    = trainingSamples;

      markDone(3);
      insertSnapshotButton();
      setStatus('Step 4 complete. Use [4b] to export covariate snapshot.');
      if (typeof onDone === 'function') onDone();
    });
  });
}


// ─────────────────────────────────────────────────────────────────
// STEP 5 — TRAIN FINAL AGBD MODEL + 3-FOLD CROSS-VALIDATION
//
// After training the final RF, a 3-fold CV provides a held-out RMSE
// that is used as the spatially constant σ for the GEDI RF in the
// inverse-variance weighted forest ensemble (Step 7).
//
// Interpreting OOB / CV RMSE for this AOI:
//   The RMSE is computed in AGBD units (Mg/ha), then converted to
//   forest carbon (kg/m²) via AGB_TO_C_KGM2 = (1+BGB_RATIO) ×
//   CARBON_FRACTION × MGHA_TO_KGM2 (= 0.061). GEDI L4A is aboveground
//   BIOMASS density, so the carbon fraction and BGB uplift are required
//   to match Sothe FC (total biomass carbon).
//   A lower RMSE = the RF predictions are closer to held-out GEDI
//   footprint values = higher weight in the ensemble.
//   Note: GEDI L4A AGBD has its own prediction error (agbd_se).
//   The CV RMSE captures RF model error *given* GEDI as truth, not
//   the absolute uncertainty in the true forest carbon stock.
//   Field sampling provides ground truth to update all estimates.
// ─────────────────────────────────────────────────────────────────
function step5_trainFinalGEDI(onDone) {
  if (!covariates)    { setStatus('Run Step 4 first - covariates not built.'); return; }
  if (!TOP_BANDS_EE)  { setStatus('Run Step 4 first - variable selection not done.'); return; }
  if (!gedi_training) { setStatus('Run Step 4 first - GEDI training data not loaded.'); return; }
  setStatus('Step 5: Training final RF (' + FINAL_N_TREES + ' trees, top ' + TOP_N_BANDS + ' bands)...');

  gedi_training.aggregate_stats('agbd').evaluate(function(s) {
    print('GEDI AGBD training distribution (Mg/ha): n=' + s.total_count +
          ' mean=' + s.mean.toFixed(2) + ' sd=' + s.total_sd.toFixed(2));
  });

  var finalRF = ee.Classifier.smileRandomForest({
    numberOfTrees: FINAL_N_TREES, minLeafPopulation: 5, bagFraction: 0.632, seed: 42
  }).setOutputMode('REGRESSION')
    .train({ features: gedi_training, classProperty: 'agbd', inputProperties: TOP_BANDS_EE });

  ee.Dictionary(finalRF.explain().get('importance')).evaluate(function(impObj) {
    var chartData = [];
    for (var b in impObj) chartData.push({ band: b, importance: impObj[b] });
    chartData.sort(function(a, bv) { return bv.importance - a.importance; });
    print(ui.Chart.feature.byFeature(
      ee.FeatureCollection(chartData.map(function(d) { return ee.Feature(null, d); })),
      'band', 'importance').setChartType('BarChart').setOptions({
        title: 'Final RF (' + FINAL_N_TREES + ' trees) - Variable Importance',
        hAxis: { title: 'Importance (Gini)' }, vAxis: { title: 'Band' },
        legend: 'none', colors: ['#1565C0'], bar: { groupWidth: '80%' }
    }));
  });

  var agbd_mgha = covariates.select(TOP_BANDS_EE).classify(finalRF).clip(aoi).rename('agbd_mgha');
  agbd_mgha     = agbd_mgha.updateMask(agbd_mgha.gte(0).and(agbd_mgha.lte(600)));
  agbd_modeled  = agbd_mgha.rename('agbd_modeled');
  // GEDI AGBD (Mg/ha) → forest carbon (kg C/m², AGB+BGB) via AGB_TO_C_KGM2
  forest_rf_pred = agbd_mgha.multiply(AGB_TO_C_KGM2).rename('forest_carbon_kgm2');

  agbd_mgha.reduceRegion({
    reducer: ee.Reducer.mean().combine(ee.Reducer.min(), null, true)
              .combine(ee.Reducer.max(), null, true).combine(ee.Reducer.stdDev(), null, true),
    geometry: aoi, scale: EXPORT_SCALE, crs: EXPORT_CRS, maxPixels: 1e11, bestEffort: true
  }).evaluate(function(r) {
    var f = function(k) { return r[k] !== undefined ? r[k].toFixed(2) : 'null'; };
    print('Forest RF: mean=' + f('agbd_mgha_mean') + ' min=' + f('agbd_mgha_min') +
          ' max=' + f('agbd_mgha_max') + ' sd=' + f('agbd_mgha_stdDev') + ' Mg/ha');
  });

  Map.addLayer(agbd_mgha, AGBD_VIS, 'AGBD - Final RF (Mg/ha)', false);
  addLayerWithStats(forest_rf_pred, 'forest_carbon_kgm2',
    { min: 0, max: 20, palette: ['#f7fcf5', '#74c476', '#00441b'] },
    'Forest Carbon - GEDI RF (kg/m2)', aoi, EXPORT_SCALE);

  Export.image.toDrive({ image: agbd_mgha.toFloat(), description: 'AGBD_RF_Top6_Extrapolation_MgHa',
    folder: EXPORT_FOLDER, region: aoi, scale: EXPORT_SCALE, crs: EXPORT_CRS, maxPixels: 1e13 });

  // ── 3-Fold Cross-Validation ────────────────────────────────────
  // Assign each training point to one of 3 folds using randomColumn.
  // For each fold: train RF on the other 2 folds, predict on the hold-out.
  // RMSE is computed server-side via sq_err → aggregate_mean → sqrt.
  // ──────────────────────────────────────────────────────────────
  var trainingWithFold = gedi_training.randomColumn('fold_rand', 42);

  var computeFoldRMSE = function(foldIdx, lo, hi) {
    var test  = trainingWithFold
      .filter(ee.Filter.gte('fold_rand', lo))
      .filter(ee.Filter.lt('fold_rand', hi));
    var train = trainingWithFold.filter(
      ee.Filter.or(ee.Filter.lt('fold_rand', lo), ee.Filter.gte('fold_rand', hi))
    );
    var rf_cv = ee.Classifier.smileRandomForest({
      numberOfTrees: FINAL_N_TREES, minLeafPopulation: 5, bagFraction: 0.632, seed: 42 + foldIdx
    }).setOutputMode('REGRESSION')
      .train({ features: train, classProperty: 'agbd', inputProperties: TOP_BANDS_EE });
    var preds = test.classify(rf_cv, 'agbd_pred');
    var sq_err = preds.map(function(f) {
      return f.set('sq_err', ee.Number(f.get('agbd_pred')).subtract(ee.Number(f.get('agbd'))).pow(2));
    });
    return ee.Number(sq_err.aggregate_mean('sq_err')).sqrt();
  };

  var rmse0_ee = computeFoldRMSE(0, 0,     1/3);
  var rmse1_ee = computeFoldRMSE(1, 1/3,   2/3);
  var rmse2_ee = computeFoldRMSE(2, 2/3,   1.0);

  ee.List([rmse0_ee, rmse1_ee, rmse2_ee]).evaluate(function(rmses, err) {
    if (err) { print('3-fold CV error: ' + err); return; }

    var r0 = rmses[0], r1 = rmses[1], r2 = rmses[2];
    var mean_rmse_mgha = (r0 + r1 + r2) / 3.0;
    // Convert AGBD RMSE (Mg/ha) → carbon σ (kg C/m², AGB+BGB), same factor as the mean
    var mean_rmse_kgm2 = mean_rmse_mgha * AGB_TO_C_KGM2;

    GEDI_RF_RMSE_MGHA   = mean_rmse_mgha;
    GEDI_RF_SIGMA_KGPM2 = mean_rmse_kgm2;

    print('');
    print('══ 3-FOLD CROSS-VALIDATION — GEDI RF ══════════════════════════');
    print('  Fold 0  RMSE: ' + r0.toFixed(2) + ' Mg/ha AGBD');
    print('  Fold 1  RMSE: ' + r1.toFixed(2) + ' Mg/ha AGBD');
    print('  Fold 2  RMSE: ' + r2.toFixed(2) + ' Mg/ha AGBD');
    print('  ──────────────────────────────────────────────────────────────');
    print('  Mean 3-fold CV RMSE:  ' + mean_rmse_mgha.toFixed(2) + ' Mg/ha AGBD');
    print('  Forest carbon σ:      ' + mean_rmse_kgm2.toFixed(3) + ' kg C/m²  (RMSE × ' + AGB_TO_C_KGM2.toFixed(3) + ', AGB+BGB)');
    print('');
    print('  Interpretation:');
    print('    The CV RMSE measures how well the RF predicts held-out GEDI');
    print('    footprint AGBD values — it is the model prediction error,');
    print('    not the uncertainty in the true forest carbon stock.');
    print('    A lower RMSE → higher weight in the Step 7 inverse-variance');
    print('    weighted ensemble relative to Sothe FC and SCANFI.');
    print('    Field measurements will update and replace all three priors.');
    print('══════════════════════════════════════════════════════════════════');
    print('');

    markDone(4);
    setStatus('Step 5 complete — CV RMSE: ' + mean_rmse_kgm2.toFixed(3) + ' kg/m². Export queued.');
    if (typeof onDone === 'function') onDone();
  });
}


// ─────────────────────────────────────────────────────────────────
// STEP 6 — SOIL SOURCE VERIFICATION
// ─────────────────────────────────────────────────────────────────
function step6_buildSoilModel(onDone) {
  if (!sothe_sc)  { setStatus('Run Step 2 first - Sothe SC not loaded.'); return; }
  if (!sg_soc_1m) { setStatus('Run Step 2 first - SoilGrids not loaded.'); return; }

  print('Step 6: Soil ensemble — Sothe et al. (2022) + SoilGrids v2.0 (0-100 cm)');
  print('  Per-pixel σ: Sothe SC unc (kg/m²)');
  print('  Derived σ:   SoilGrids — standardized(Sothe unc) + standardized(inter-product SD)');
  print('               scaled to Sothe mean σ. Computed in Step 7.');

  var soilVis = { min: 0, max: 30, palette: ['#fff7bc', '#fe9929', '#993404'] };
  addLayerWithStats(sothe_sc,  'sothe_sc',  soilVis, 'Sothe SC - soil source 1 (kg/m2)', aoi, 250);
  addLayerWithStats(sg_soc_1m, 'sg_soc_1m', soilVis, 'SoilGrids 0-100 cm - soil source 2 (kg/m2)', aoi, 250);

  var diff = sothe_sc.subtract(sg_soc_1m).rename('sothe_minus_sg');
  Map.addLayer(diff, { min: -15, max: 15, palette: ['#d73027','#ffffbf','#1a9850'] },
    'Soil Source Diff: Sothe - SoilGrids (kg/m2)', false);
  printStats(diff, 'sothe_minus_sg', aoi, 250, 'Soil Source Diff: Sothe - SoilGrids (kg/m2)');

  markDone(5);
  setStatus('Step 6 complete - soil sources verified. Run Step 7 to build ensemble.');
  if (typeof onDone === 'function') onDone();
}


// ─────────────────────────────────────────────────────────────────
// STEP 7 — INVERSE-VARIANCE WEIGHTED ENSEMBLE MODELS
//
// Replaces equal-weight ensemble mean with precision-weighted means.
//
// FOREST — 3 members:
//   GEDI RF:  σ = GEDI_RF_SIGMA_KGPM2 (constant, from 3-fold CV)
//   Sothe FC: σ = sothe_fc_unc per-pixel (kg/m²)
//   SCANFI:   σ = SCANFI_SIGMA_FC = 1.935 kg/m² (constant, Guindon 2024)
//   w_i(p) = 1/σ_i(p)²
//   forest_wmean(p) = Σ[w_i(p) × x_i(p)] / Σ[w_i(p)]
//   forest_wmean_sigma(p) = 1 / sqrt(Σ[w_i(p)])
//
// SOIL — 2 members:
//   Sothe SC: σ = sothe_sc_unc per-pixel (kg/m²)
//   SoilGrids: σ derived (no published per-pixel layer in GEE yet):
//     A = sothe_sc_unc / mean(sothe_sc_unc)  — dimensionless, mean=1
//     B = soil_ens_sd  / mean(soil_ens_sd)   — dimensionless, mean=1
//     combined_index   = (A + B) / 2         — equal-weight average
//     sg_sigma = combined_index × mean(sothe_sc_unc)
//     Rationale: SoilGrids is most uncertain where Sothe itself is
//     uncertain AND where the two products disagree. Standardising
//     before averaging prevents either component from dominating by
//     magnitude. Scaling back to Sothe mean σ anchors the result to
//     a physically meaningful baseline.
//
//   NOTE: The two normalisation means (sothe_sc_unc mean, soil_ens_sd
//   mean) are resolved to plain JS numbers via evaluate() before any
//   image operations run. Using lazy ee.Number(dict.get()) inside
//   image.divide() causes a 'key not found' error in GEE's graph
//   evaluator when the dict hasn't been computed yet.
//
// Equal-weight ensemble means are retained globally for comparison.
// forest_nf_mask updated to use forest_wmean.
// total_ecosystem_c_w = forest_wmean + soil_wmean.
// ─────────────────────────────────────────────────────────────────
function step7_buildEnsemble(onDone) {
  if (!forest_rf_pred)      { setStatus('Run Step 5 first - forest RF missing.'); return; }
  if (!GEDI_RF_SIGMA_KGPM2) { setStatus('Run Step 5 first - 3-fold CV RMSE not computed.'); return; }
  if (!sothe_sc)            { setStatus('Run Step 2 first - Sothe SC missing.'); return; }
  if (!sg_soc_1m)           { setStatus('Run Step 2 first - SoilGrids missing.'); return; }
  if (!sothe_fc)            { setStatus('Run Step 2 first - Sothe FC missing.'); return; }
  if (!sothe_fc_unc)        { setStatus('Run Step 2 first - Sothe FC unc missing.'); return; }
  if (!sothe_sc_unc)        { setStatus('Run Step 2 first - Sothe SC unc missing.'); return; }
  if (!scanfi_fc)           { setStatus('Run Step 3 first - SCANFI FC missing.'); return; }
  if (!sbfi_agb_raster)     { setStatus('Run Step 4 first - SBFI not rasterised.'); return; }

  setStatus('Step 7: Building equal-weight reference ensemble...');

  // ── Equal-weight forest ensemble (kept for comparison) ─────────
  // v4.6: uses the SAME 3 members as the weighted ensemble (GEDI RF, Sothe FC,
  // SCANFI) so "weighted vs equal-weight" is a like-for-like comparison.
  // (was GEDI RF + Sothe FC + SBFI, which had no σ for the weighted version.)
  var forest_eq_stack = forest_rf_pred.rename('gedi_rf')
    .addBands(sothe_fc.rename('sothe_fc'))
    .addBands(scanfi_fc.rename('scanfi_fc'));
  var f_count   = forest_eq_stack.reduce(ee.Reducer.count()).rename('f_model_count');
  var f_mask_eq = f_count.gte(2);
  forest_ens_mean = forest_eq_stack.reduce(ee.Reducer.mean()).updateMask(f_mask_eq).rename('forest_ens_mean');
  forest_ens_sd   = forest_eq_stack.reduce(ee.Reducer.stdDev()).updateMask(f_mask_eq).rename('forest_ens_sd');

  // ── Equal-weight soil ensemble (kept for comparison + sg_sigma derivation) ──
  var soil_stack = sothe_sc.rename('sothe').addBands(sg_soc_1m.rename('soilgrids'));
  var s_count    = soil_stack.reduce(ee.Reducer.count()).rename('s_model_count');
  var s_mask_eq  = s_count.gte(2);
  soil_ens_mean = soil_stack.reduce(ee.Reducer.mean()).updateMask(s_mask_eq).rename('soil_ens_mean');
  soil_ens_sd   = soil_stack.reduce(ee.Reducer.stdDev()).updateMask(s_mask_eq).rename('soil_ens_sd');

  // ── Legacy RSS uncertainty (retained for Neyman allocation) ────
  var gedi_se_pct = gedi_agbd_se.divide(forest_ens_mean.abs().add(0.001));
  var raw_f_sd    = forest_ens_sd.where(gedi_agbd_se.mask(), forest_ens_sd.add(gedi_se_pct).divide(2));
  forest_uncertainty = raw_f_sd.pow(2).add(sothe_fc_unc.pow(2))
    .sqrt().updateMask(f_mask_eq).rename('forest_uncertainty_rss');
  soil_uncertainty = sothe_sc_unc.pow(2).add(soil_ens_sd.pow(2))
    .sqrt().updateMask(s_mask_eq).rename('soil_uncertainty_rss');

  // ── Inverse-variance weighted forest ensemble ──────────────────
  // v4.6: members that are absent at a pixel contribute zero weight (rather
  // than nulling the whole sum). Each weighted term and its weight are set to
  // 0 where the member is masked, so Σ[w·x] / Σ[w] is taken over only the
  // members actually present. The result is then masked to require ≥2 members
  // (f_mask_eq), so the "≥2 of 3" rule is genuinely honoured.
  setStatus('Step 7: Building inverse-variance weighted ensembles...');

  var sigma_gedi_img   = ee.Image.constant(GEDI_RF_SIGMA_KGPM2).rename('sigma_gedi');
  var sigma_sothe_f    = sothe_fc_unc.rename('sigma_sothe_f');
  var sigma_scanfi_img = ee.Image.constant(SCANFI_SIGMA_FC).rename('sigma_scanfi');

  // Per-member presence masks
  var m_gedi   = forest_rf_pred.mask();
  var m_sothe  = sothe_fc.mask();
  var m_scanfi = scanfi_fc.mask();

  // Effective weights (0 where member absent)
  var we_gedi   = sigma_gedi_img.pow(2).pow(-1).multiply(m_gedi).unmask(0);
  var we_sothe  = sigma_sothe_f.pow(2).pow(-1).multiply(m_sothe).unmask(0);
  var we_scanfi = sigma_scanfi_img.pow(2).pow(-1).multiply(m_scanfi).unmask(0);

  // Weighted value terms (0 where member absent)
  var wx_gedi   = forest_rf_pred.unmask(0).multiply(we_gedi);
  var wx_sothe  = sothe_fc.unmask(0).multiply(we_sothe);
  var wx_scanfi = scanfi_fc.unmask(0).multiply(we_scanfi);

  var w_sum_f = we_gedi.add(we_sothe).add(we_scanfi);

  forest_wmean = wx_gedi.add(wx_sothe).add(wx_scanfi)
    .divide(w_sum_f)
    .updateMask(f_mask_eq)
    .rename('forest_wmean');

  forest_wmean_sigma = w_sum_f.pow(-0.5).updateMask(f_mask_eq).rename('forest_wmean_sigma');
  forest_nf_mask     = forest_wmean.gt(FOREST_THRESHOLD);

  // Data-source map: number of forest products present per pixel (1–3)
  dsm_forest = f_count.clip(aoi).rename('forest_data_source');
  dsm_soil = ee.Image(1)
    .where(sothe_sc.mask().not(),   ee.Image(2))
    .where(sg_soc_1m.mask().not(), ee.Image(3))
    .where(sothe_sc.mask().not().and(sg_soc_1m.mask().not()), ee.Image(4))
    .clip(aoi).rename('soil_data_source');
  total_ecosystem_c = forest_ens_mean.add(soil_ens_mean).rename('total_ec_equal_weight');

  // Forest pixel count + console
  forest_wmean.reduceRegion({ reducer: ee.Reducer.count(), geometry: aoi,
    scale: EXPORT_SCALE, crs: EXPORT_CRS, maxPixels: 1e11, bestEffort: true })
    .evaluate(function(r) { print('Valid forest weighted-mean pixels:', r['forest_wmean'] || 0); });

  print('');
  print('══ INVERSE-VARIANCE WEIGHTED ENSEMBLE — Forest ═══════════════');
  print('  σ per product:');
  print('    GEDI RF:  [constant]  ' + GEDI_RF_SIGMA_KGPM2.toFixed(3) + ' kg/m²  (3-fold CV RMSE × ' + AGB_TO_C_KGM2.toFixed(3) + ', AGB+BGB)');
  print('    Sothe FC: [per-pixel] see map layer "Sothe FC Uncertainty"');
  print('    SCANFI:   [constant]  ' + SCANFI_SIGMA_FC.toFixed(3) + ' kg/m²  (Guindon 2024 RMSE × ' + AGB_TO_C_KGM2.toFixed(3) + ')');
  print('  Formula: w_i(p) = 1/σ_i(p)²  |  wmean(p) = Σ[w_i·x_i] / Σ[w_i]');
  print('══════════════════════════════════════════════════════════════════');
  print('');

  printStats(forest_wmean,       'forest_wmean',       aoi, EXPORT_SCALE, 'Forest Weighted Mean (kg/m2)');
  printStats(forest_wmean_sigma, 'forest_wmean_sigma', aoi, EXPORT_SCALE, 'Forest Weighted σ (kg/m2)');
  printStats(forest_ens_mean,    'forest_ens_mean',    aoi, EXPORT_SCALE, 'Forest Equal-weight Mean [comparison] (kg/m2)');

  var carbonVis_f = { min: 0, max: 20, palette: ['#f7f7f7', '#2ca25f', '#006837'] };
  var sigmaVis_f  = { min: 0, max: 8,  palette: ['#2166ac', '#f7f7f7', '#d6604d'] };
  Map.addLayer(forest_wmean,       carbonVis_f, 'Forest Weighted Mean (kg/m2)',                  true);
  Map.addLayer(forest_wmean_sigma, sigmaVis_f,  'Forest Weighted σ (kg/m2)',                     false);
  Map.addLayer(forest_ens_mean,    carbonVis_f, 'Forest Equal-weight Mean [comparison] (kg/m2)', false);
  Map.addLayer(forest_ens_sd,      sigmaVis_f,  'Forest Model Spread SD (kg/m2)',                false);
  Map.addLayer(sothe_fc_unc,       sigmaVis_f,  'Sothe FC Uncertainty / per-pixel (kg/m2)',      false);
  Map.addLayer(
    forest_wmean.subtract(forest_ens_mean).rename('forest_w_minus_eq'),
    { min: -5, max: 5, palette: ['#d73027', '#ffffbf', '#1a9850'] },
    'Forest: Weighted − Equal-weight Mean (kg/m2)', false);
  Map.addLayer(forest_nf_mask, { min: 0, max: 1, palette: ['#ffffff', '#006837'] },
    'Forest / Non-Forest Mask (weighted mean)', false);

  // ── Soil: resolve normalisation means client-side before image ops ──
  // Two sequential evaluate() calls are used instead of ee.List([dict.get(), ...])
  // because GEE does not reliably serialise ee.Object results from .get() inside
  // an ee.List for evaluate(). Evaluating each reduceRegion dict independently
  // avoids the 'key not found' error and makes failures easier to diagnose.
  sothe_sc_unc.reduceRegion({
    reducer: ee.Reducer.mean(), geometry: aoi, scale: 250,
    crs: EXPORT_CRS, maxPixels: 1e11, bestEffort: true
  }).evaluate(function(r1, e1) {

    if (e1 || !r1) {
      print('ERROR reducing sothe_sc_unc: ' + (e1 || 'null result'));
      setStatus('Step 7 soil error — see console.'); return;
    }
    // Diagnostic: print all keys so band-name mismatches are immediately visible
    print('sothe_sc_unc reduceRegion keys:', Object.keys(r1));

    var sothe_unc_mean_js = r1['sothe_sc_unc'];
    if (sothe_unc_mean_js === undefined || sothe_unc_mean_js === null) {
      print('ERROR: sothe_sc_unc not in result. Keys returned:', Object.keys(r1));
      setStatus('Step 7 soil error — unexpected band name. See console.'); return;
    }

    soil_ens_sd.reduceRegion({
      reducer: ee.Reducer.mean(), geometry: aoi, scale: 250,
      crs: EXPORT_CRS, maxPixels: 1e11, bestEffort: true
    }).evaluate(function(r2, e2) {

      if (e2 || !r2) {
        print('ERROR reducing soil_ens_sd: ' + (e2 || 'null result'));
        setStatus('Step 7 soil error — see console.'); return;
      }
      print('soil_ens_sd reduceRegion keys:', Object.keys(r2));

      var soil_sd_mean_js = r2['soil_ens_sd'];
      if (soil_sd_mean_js === undefined || soil_sd_mean_js === null) {
        print('ERROR: soil_ens_sd not in result. Keys returned:', Object.keys(r2));
        setStatus('Step 7 soil error — unexpected band name. See console.'); return;
      }

    print('══ INVERSE-VARIANCE WEIGHTED ENSEMBLE — Soil ═════════════════');
    print('  Normalisation means (client-side resolved):');
    print('    Sothe SC unc mean:      ' + sothe_unc_mean_js.toFixed(3) + ' kg/m²  [key: sothe_sc_unc]');
    print('    Inter-product SD mean:  ' + soil_sd_mean_js.toFixed(3)   + ' kg/m²  [key: soil_ens_sd]');
    print('  SoilGrids σ derivation:');
    print('    A = Sothe SC unc / ' + sothe_unc_mean_js.toFixed(3) + '  [dimensionless, mean=1]');
    print('    B = inter-product SD / ' + soil_sd_mean_js.toFixed(3) + '  [dimensionless, mean=1]');
    print('    sg_sigma = ((A + B) / 2) × ' + sothe_unc_mean_js.toFixed(3) + '  [kg/m²]');
    print('    Result: SoilGrids σ is highest where Sothe is uncertain');
    print('    AND where the two products disagree most.');
    print('══════════════════════════════════════════════════════════════════');
    print('');

    // Derived SoilGrids σ using plain JS scalar divisors
    var A = sothe_sc_unc.divide(sothe_unc_mean_js);
    var B = soil_ens_sd.divide(soil_sd_mean_js);
    sg_sigma = A.add(B).divide(2).multiply(sothe_unc_mean_js)
      .rename('sg_sigma').updateMask(s_mask_eq);

    // Inverse-variance soil ensemble
    var w_sothe_s = sothe_sc_unc.pow(2).pow(-1);
    var w_sg_s    = sg_sigma.pow(2).pow(-1);
    var w_sum_s   = w_sothe_s.add(w_sg_s).updateMask(s_mask_eq);

    soil_wmean = sothe_sc.multiply(w_sothe_s)
      .add(sg_soc_1m.multiply(w_sg_s))
      .divide(w_sum_s)
      .rename('soil_wmean');

    soil_wmean_sigma    = w_sum_s.pow(-0.5).rename('soil_wmean_sigma');
    total_ecosystem_c_w = forest_wmean.add(soil_wmean).rename('total_ec_weighted');

    // Pixel count validation
    soil_wmean.reduceRegion({ reducer: ee.Reducer.count(), geometry: aoi,
      scale: 250, crs: EXPORT_CRS, maxPixels: 1e11, bestEffort: true })
      .evaluate(function(r) { print('Valid soil weighted-mean pixels:', r['soil_wmean'] || 0); });

    // Console stats — soil
    printStats(sg_sigma,            'sg_sigma',           aoi, 250,          'SoilGrids Derived σ (kg/m2)');
    printStats(soil_wmean,          'soil_wmean',         aoi, 250,          'Soil Weighted Mean (kg/m2)');
    printStats(soil_wmean_sigma,    'soil_wmean_sigma',   aoi, 250,          'Soil Weighted σ (kg/m2)');
    printStats(soil_ens_mean,       'soil_ens_mean',      aoi, 250,          'Soil Equal-weight Mean [comparison] (kg/m2)');
    printStats(total_ecosystem_c_w, 'total_ec_weighted',  aoi, EXPORT_SCALE, 'Total Ecosystem Carbon - Weighted (kg/m2)');

    // Map layers — soil and total
    var soilCarVis = { min: 15, max: 100, palette: ['#fff7bc', '#fe9929', '#993404'] };
    var sigmaVis_s = { min: 0,  max: 50,  palette: ['#2166ac', '#f7f7f7', '#d6604d'] };
    var totalVis   = { min: 0,  max: 60,  palette: ['#f7f7f7', '#2ca25f', '#006837'] };

    Map.addLayer(soil_wmean,       soilCarVis, 'Soil Weighted Mean (kg/m2)',                  true);
    Map.addLayer(soil_wmean_sigma, sigmaVis_s, 'Soil Weighted σ (kg/m2)',                     false);
    Map.addLayer(soil_ens_mean,    soilCarVis, 'Soil Equal-weight Mean [comparison] (kg/m2)', false);
    Map.addLayer(sg_sigma,         sigmaVis_s, 'SoilGrids Derived σ (kg/m2)',                 false);
    Map.addLayer(sothe_sc_unc,     sigmaVis_s, 'Sothe SC Uncertainty / per-pixel (kg/m2)',    false);
    Map.addLayer(soil_ens_sd,      sigmaVis_s, 'Soil Inter-product Disagreement SD (kg/m2)',  false);
    Map.addLayer(
      soil_wmean.subtract(soil_ens_mean).rename('soil_w_minus_eq'),
      { min: -5, max: 5, palette: ['#d73027', '#ffffbf', '#1a9850'] },
      'Soil: Weighted − Equal-weight Mean (kg/m2)', false);
    Map.addLayer(total_ecosystem_c_w, totalVis,
      'Total Ecosystem Carbon - Weighted (kg/m2)',     false);
    Map.addLayer(total_ecosystem_c,   totalVis,
      'Total Ecosystem Carbon - Equal-weight (kg/m2)', false);
    Map.addLayer(dsm_forest,
      { min: 1, max: 3, palette: ['#feb24c', '#74c476', '#006837'] },
      'Forest Data Source (# products present: 1–3)', false);
    Map.addLayer(dsm_soil,
      { min: 1, max: 4, palette: ['#2ca25f', '#feb24c', '#de2d26', '#969696'] },
      'Soil Data Source (1=Both 2=SG only 3=Sothe only 4=Neither)', false);

    markDone(6);
    setStatus('Step 7 complete - inverse-variance weighted ensembles built.');
    if (typeof onDone === 'function') onDone();

    }); // end soil_ens_sd evaluate()
  }); // end sothe_sc_unc evaluate()
}


// ─────────────────────────────────────────────────────────────────
// STEP 8 — NEYMAN-ALLOCATION STRATIFIED SAMPLING
// Uses forest_wmean and soil_wmean (weighted best estimates) for
// the Neyman allocation and power analysis, replacing equal-weight.
// RSS uncertainty images retained as stratification signal since
// they represent the full spread across all error sources.
// ─────────────────────────────────────────────────────────────────
function step8_generateSampling(onDone) {
  if (!forest_uncertainty) { setStatus('Run Step 7 first - forest uncertainty missing.'); return; }
  if (!soil_uncertainty)   { setStatus('Run Step 7 first - soil uncertainty missing.');   return; }
  if (!forest_nf_mask)     { setStatus('Run Step 7 first - forest mask missing.');        return; }
  if (!forest_wmean)       { setStatus('Run Step 7 first - forest weighted mean missing.'); return; }
  if (!soil_wmean)         { setStatus('Run Step 7 first - soil weighted mean missing.'); return; }
  setStatus('Step 8: Loading ESA WorldCover land-cover strata...');

  // Fire onDone only after BOTH forest and soil sampling have finished
  var _poolsDone = 0;
  var _afterPool = function() {
    _poolsDone += 1;
    if (_poolsDone === 2) {
      markDone(7);
      setStatus('Step 8 complete - Neyman-allocated sampling points generated.');
      if (typeof onDone === 'function') onDone();
    }
  };

  var worldcover = ee.ImageCollection('ESA/WorldCover/v200').first()
    .reproject({ crs: EXPORT_CRS, scale: EXPORT_SCALE }).clip(aoi);

  var lc_broad = worldcover.remap(
    [10,  20,  30,  40,  50,  60,  70,  80,  90,  95, 100],
    [ 1,   2,   2,   4,   5,   5,   5,   0,   3,   3,   2]
  ).rename('lc_class').toInt();
  lc_broad = lc_broad.updateMask(lc_broad.gt(0));

  Map.addLayer(lc_broad, { min: 1, max: 5,
    palette: ['#006837','#addd8e','#41b6c4','#fec44f','#bdbdbd'] },
    'ESA WorldCover Broad LC (1=Forest 2=Shrub 3=Wetland 4=Crop 5=Other)', false);

  // Extra bands carry weighted means + σ to sampled points
  var forestExtras = forest_wmean.unmask(0)
    .addBands(forest_wmean_sigma.unmask(0))
    .addBands(forest_ens_mean.unmask(0))
    .addBands(forest_ens_sd.unmask(0))
    .addBands(forest_uncertainty.unmask(0))
    .addBands(dsm_forest)
    .addBands(lc_broad);

  var soilExtras = soil_wmean.unmask(0)
    .addBands(soil_wmean_sigma.unmask(0))
    .addBands(soil_ens_mean.unmask(0))
    .addBands(soil_ens_sd.unmask(0))
    .addBands(soil_uncertainty.unmask(0))
    .addBands(sg_sigma.unmask(0))
    .addBands(dsm_soil)
    .addBands(lc_broad);

  var forest_unc_masked = forest_uncertainty.updateMask(forest_nf_mask);
  var lc_forest_masked  = lc_broad.updateMask(forest_nf_mask);

  _neymanSample(forest_unc_masked, 'forest_uncertainty_rss', lc_forest_masked, 'lc_class',
    forestExtras, N_FOREST_SAMPLES, 'Forest',
    function(pts) {
      forest_sampling_pts = pts;
      pts.size().evaluate(function(n) {
        Map.addLayer(pts, { color: '1565c0' }, 'Forest Sampling Points - Neyman (' + n + ')', true);
        print('Forest Neyman sampling complete: ' + n + ' points.');
        _afterPool();
      });
    }
  );

  _neymanSample(soil_uncertainty, 'soil_uncertainty_rss', lc_broad, 'lc_class',
    soilExtras, N_SOIL_SAMPLES, 'Soil',
    function(pts) {
      soil_sampling_pts = pts;
      pts.size().evaluate(function(n) {
        Map.addLayer(pts, { color: 'e65100' }, 'Soil Sampling Points - Neyman (' + n + ')', true);
        print('Soil Neyman sampling complete: ' + n + ' points.');
        _afterPool();
      });
    }
  );

  _powerAnalysis();
}

// ─────────────────────────────────────────────────────────────────
// _neymanSample — Neyman optimal allocation (unchanged from v4.4)
// ─────────────────────────────────────────────────────────────────
function _neymanSample(unc_img, unc_band, lc_img, lc_band, extra_bands, n_total, pool_label, callback) {
  var unc = unc_img.select(unc_band);

  unc.reduceRegion({ reducer: ee.Reducer.percentile([25, 50, 75]),
    geometry: aoi, scale: EXPORT_SCALE, crs: EXPORT_CRS,
    maxPixels: 1e11, bestEffort: true
  }).evaluate(function(pct, pctErr) {
    if (pctErr) { print('ERROR computing percentiles (' + pool_label + '): ' + pctErr); return; }

    var p25 = pct[unc_band + '_p25'] || 1;
    var p50 = pct[unc_band + '_p50'] || 3;
    var p75 = pct[unc_band + '_p75'] || 6;

    print(pool_label + ' uncertainty quartile breaks:');
    print('  p25=' + p25.toFixed(4) + ' | p50=' + p50.toFixed(4) + ' | p75=' + p75.toFixed(4));

    var binSigmas = [p25/2, (p25+p50)/2, (p50+p75)/2, p75*1.5];
    print('  Bin σ_h midpoints: [' + binSigmas.map(function(v){ return v.toFixed(4); }).join(', ') + ']');

    var unc_unm = unc.unmask(0);
    var unc_bin = ee.Image(N_UNC_BINS - 1)
      .where(unc_unm.lt(p75), N_UNC_BINS - 2)
      .where(unc_unm.lt(p50), N_UNC_BINS - 3)
      .where(unc_unm.lt(p25), 0)
      .toInt().rename('unc_bin');

    var stratum = lc_img.select(lc_band).multiply(N_UNC_BINS).add(unc_bin).toInt().rename('stratum');

    Map.addLayer(unc_bin.updateMask(lc_img.select(lc_band).mask()),
      { min: 0, max: N_UNC_BINS-1, palette: ['#2c7bb6','#abd9e9','#fdae61','#d7191c'] },
      pool_label + ' Uncertainty Bins (0=low 3=high)', false);
    Map.addLayer(stratum,
      { min: N_UNC_BINS, max: 5*N_UNC_BINS + N_UNC_BINS - 1,
        palette: ['#d4e6f1','#85c1e9','#2e86c1','#1a5276',
                  '#d5f5e3','#76d7c4','#1abc9c','#0e6655',
                  '#e8daef','#bb8fce','#8e44ad','#6c3483',
                  '#fef9e7','#f8c471','#e67e22','#ca6f1e',
                  '#f2f3f4','#aab7b8','#717d7e','#424949'] },
      pool_label + ' Composite Strata (LC × unc_bin)', false);

    stratum.reduceRegion({ reducer: ee.Reducer.frequencyHistogram(),
      geometry: aoi, scale: EXPORT_SCALE, crs: EXPORT_CRS,
      maxPixels: 1e11, bestEffort: true
    }).get('stratum').evaluate(function(histObj, histErr) {
      if (histErr) { print('ERROR stratum histogram (' + pool_label + '): ' + histErr); return; }
      if (!histObj || Object.keys(histObj).length === 0) {
        setStatus(pool_label + ': empty histogram - check AOI coverage.'); return;
      }

      var keys = Object.keys(histObj), totalWeight = 0, strata = {};
      keys.forEach(function(h) {
        var hInt = parseInt(h), lc_code = Math.floor(hInt/N_UNC_BINS), unc_b = hInt % N_UNC_BINS;
        if (lc_code < 1) return;
        var N_h = histObj[h] || 0, sigma_h = binSigmas[unc_b] || 0, w_h = N_h * sigma_h;
        strata[h] = { N: N_h, sigma: sigma_h, w: w_h, lc: lc_code, bin: unc_b };
        totalWeight += w_h;
      });

      var classValues = [], classPoints = [];
      var lcNames = ['?','Forest','Shrub/Grass','Wetland','Cropland','Other'];
      print('NEYMAN ALLOCATION - ' + pool_label);
      print('  Total Neyman weight Σ(N_h × σ_h) = ' + totalWeight.toFixed(2));

      // ── Largest-remainder allocation → sums to EXACTLY n_total (v4.6) ──
      // Eligible strata = occupied (N_h>0) with positive Neyman weight.
      // Each gets floor(ideal); the leftover points go to the largest
      // fractional remainders. This guarantees Σ n_h = n_total, unlike the
      // previous per-stratum max(MIN, round(...)) which could overshoot.
      var elig = [];
      keys.forEach(function(h) {
        var s = strata[h];
        if (s && s.w > 0 && s.N > 0) elig.push({ h: parseInt(h), s: s, base: 0, rem: 0 });
      });

      if (elig.length === 0 || totalWeight <= 0) {
        setStatus(pool_label + ': no valid strata.'); return;
      }

      var assigned = 0;
      elig.forEach(function(e) {
        var ideal = n_total * e.s.w / totalWeight;
        e.base = Math.floor(ideal);
        e.rem  = ideal - e.base;
        assigned += e.base;
      });
      var remaining = n_total - assigned;
      // Hand out remaining points to the largest fractional remainders, looping
      // if n_total exceeds the stratum count (so totals always reconcile).
      elig.sort(function(a, b) { return b.rem - a.rem; });
      var gi = 0;
      while (remaining > 0) { elig[gi % elig.length].base += 1; remaining--; gi++; }

      // Emit in stratum order; keep only strata that received ≥1 point
      var nZero = 0, nBelowFloor = 0;
      elig.sort(function(a, b) { return a.h - b.h; });
      print('  Strata (largest-remainder, exact total):');
      elig.forEach(function(e) {
        var s = e.s;
        if (e.base === 0) { nZero++; return; }
        if (e.base < MIN_PTS_PER_STRATUM) nBelowFloor++;
        print('    h=' + e.h + ' (' + (lcNames[s.lc]||'LC'+s.lc) + ', bin=' + s.bin + ')' +
              '  N_h=' + s.N.toFixed(0) + '  σ_h=' + s.sigma.toFixed(4) + '  n_h=' + e.base);
        classValues.push(e.h);
        classPoints.push(e.base);
      });
      print('  Total allocated = ' + classPoints.reduce(function(a,b){return a+b;},0) +
            '  (target = ' + n_total + ', EXACT)');
      if (nZero > 0)
        print('  NOTE: ' + nZero + ' occupied stratum/strata received 0 points — ' +
              'n_total is small relative to the number of strata. Increase ' +
              'N or reduce N_UNC_BINS to spread coverage.');
      if (nBelowFloor > 0)
        print('  NOTE: ' + nBelowFloor + ' stratum/strata below the advisory floor of ' +
              MIN_PTS_PER_STRATUM + ' points (variance there will be poorly estimated).');

      if (classValues.length === 0) { setStatus(pool_label + ': no valid strata.'); return; }

      var pts = extra_bands.addBands(stratum)
        .stratifiedSample({ numPoints: 0, classBand: 'stratum', region: aoi,
          scale: EXPORT_SCALE, classValues: classValues, classPoints: classPoints,
          geometries: true, seed: 42, tileScale: 4, dropNulls: true })
        .map(function(f) { return f.set('pool', pool_label); });

      callback(pts);
    });
  });
}

// ─────────────────────────────────────────────────────────────────
// _powerAnalysis — uses weighted means (spatial SD from wmean)
// ─────────────────────────────────────────────────────────────────
function _powerAnalysis() {
  if (!forest_wmean || !soil_wmean || !forest_nf_mask) {
    print('Power analysis skipped — run Steps 5-7 first.'); return;
  }

  var f_stats = forest_wmean.updateMask(forest_nf_mask).reduceRegion({
    reducer: ee.Reducer.mean().combine(ee.Reducer.stdDev(), null, true),
    geometry: aoi, scale: EXPORT_SCALE, crs: EXPORT_CRS, maxPixels: 1e11, bestEffort: true
  });
  var s_stats = soil_wmean.reduceRegion({
    reducer: ee.Reducer.mean().combine(ee.Reducer.stdDev(), null, true),
    geometry: aoi, scale: 250, crs: EXPORT_CRS, maxPixels: 1e11, bestEffort: true
  });

  ee.List([f_stats.get('forest_wmean_mean'), f_stats.get('forest_wmean_stdDev'),
           s_stats.get('soil_wmean_mean'),   s_stats.get('soil_wmean_stdDev')
  ]).evaluate(function(vals, err) {
    if (err || !vals || vals.some(function(v){return v===null;})) {
      print('Power analysis error: ' + (err||'null values')); return;
    }
    var f_mean=vals[0], f_sd=vals[1], s_mean=vals[2], s_sd=vals[3];
    var f_cv = f_sd / Math.max(Math.abs(f_mean), 0.001);
    var s_cv = s_sd / Math.max(Math.abs(s_mean), 0.001);
    var t90 = 1.645, MOE_TARGET = 0.20, N_CHART = 100;
    var f_nmin = Math.ceil(Math.pow(t90*f_cv/MOE_TARGET, 2));
    var s_nmin = Math.ceil(Math.pow(t90*s_cv/MOE_TARGET, 2));

    print('');
    print('════ STATISTICAL POWER ANALYSIS (Step 8) — Weighted Means ══════');
    print('Estimand:  Mean carbon density (kg/m²)  |  Basis: weighted ensemble');
    print('CI: 90%  |  MOE benchmark: ±20% of mean  |  t = 1.645');
    print('');
    print('FOREST (forested pixels only, weighted mean):');
    print('  μ = ' + f_mean.toFixed(3) + ' kg/m²  |  σ = ' + f_sd.toFixed(3) +
          ' kg/m²  |  CV = ' + (f_cv*100).toFixed(1) + '%');
    print('  n_min for ±20% MOE (SRS):  ' + f_nmin);
    print('  Configured n = ' + N_FOREST_SAMPLES + '  →  expected MOE = ±' +
          (100*t90*f_cv/Math.sqrt(N_FOREST_SAMPLES)).toFixed(1) + '%' +
          ((100*t90*f_cv/Math.sqrt(N_FOREST_SAMPLES)) <= 20 ? '  ✓' : '  ⚠'));
    print('  RECOMMENDATION: ' + (N_FOREST_SAMPLES >= f_nmin
          ? 'configured n meets the ±20% target.'
          : 'increase N_FOREST_SAMPLES to ≥ ' + f_nmin + ' to meet ±20% (SRS bound).'));
    print('');
    print('SOIL (full AOI, weighted mean):');
    print('  μ = ' + s_mean.toFixed(3) + ' kg/m²  |  σ = ' + s_sd.toFixed(3) +
          ' kg/m²  |  CV = ' + (s_cv*100).toFixed(1) + '%');
    print('  n_min for ±20% MOE (SRS):  ' + s_nmin);
    print('  Configured n = ' + N_SOIL_SAMPLES + '  →  expected MOE = ±' +
          (100*t90*s_cv/Math.sqrt(N_SOIL_SAMPLES)).toFixed(1) + '%' +
          ((100*t90*s_cv/Math.sqrt(N_SOIL_SAMPLES)) <= 20 ? '  ✓' : '  ⚠'));
    print('  RECOMMENDATION: ' + (N_SOIL_SAMPLES >= s_nmin
          ? 'configured n meets the ±20% target.'
          : 'increase N_SOIL_SAMPLES to ≥ ' + s_nmin + ' to meet ±20% (SRS bound).'));
    print('');
    print('  CAVEAT: n_min uses the map\'s pixel-to-pixel spatial SD as a proxy');
    print('  for field-plot variance and assumes simple random sampling. Neyman');
    print('  stratification typically beats this bound, so these n are');
    print('  conservative. Field data will give the definitive MOE.');
    print('═════════════════════════════════════════════════════════════════');

    var features = [];
    for (var n = 1; n <= N_CHART; n++) {
      features.push(ee.Feature(null, {
        'n': n,
        'Forest MOE (%)': 100*t90*f_cv/Math.sqrt(n),
        'Soil MOE (%)':   100*t90*s_cv/Math.sqrt(n),
        'Target ±20%':    20
      }));
    }
    print(ui.Chart.feature.byFeature(ee.FeatureCollection(features), 'n',
      ['Forest MOE (%)', 'Soil MOE (%)', 'Target ±20%'])
      .setChartType('LineChart').setOptions({
        title: 'Sampling Power: MOE vs. n (90% CI, SRS bound) — Weighted Ensemble',
        hAxis: { title: 'Field samples collected (n)', viewWindow: {min:1, max:N_CHART} },
        vAxis: { title: 'Relative MOE (%)', viewWindow: {min:0, max:100} },
        series: {
          0: { color: '#1565c0', lineWidth: 2.5, pointSize: 0 },
          1: { color: '#e65100', lineWidth: 2.5, pointSize: 0 },
          2: { color: '#2e7d32', lineWidth: 1.5, lineDashStyle: [8,4], pointSize: 0 }
        }, legend: { position: 'bottom' }
      }));
  });
}


// ─────────────────────────────────────────────────────────────────
// STEP 9 — REPORTS AND EXPORTS
//
// v4.5 additions:
//   Model performance table: per-product AOI stats + σ type/value.
//   Total carbon stocks in tonnes (Mg C) for forest, soil, total.
//   Weighted means used as primary estimates; equal-weight retained
//   for comparison.
// ─────────────────────────────────────────────────────────────────
function step9_reportsAndExports(onDone) {
  if (!forest_wmean)       { setStatus('Run Step 7 first.'); return; }
  if (!soil_wmean)         { setStatus('Run Step 7 first.'); return; }
  if (!total_ecosystem_c_w){ setStatus('Run Step 7 first.'); return; }
  setStatus('Step 9: Computing AOI-wide summary statistics...');

  // ── Per-product AOI stats ──────────────────────────────────────
  print('');
  print('══ MODEL PERFORMANCE — WITHIN AOI ════════════════════════════════');
  print('');
  print('FOREST CARBON PRODUCTS  (σ = per-product uncertainty used in ensemble weighting)');
  printStats(forest_rf_pred,   'forest_carbon_kgm2', aoi, EXPORT_SCALE, 'GEDI RF (this study, kg/m2)');
  printStats(sothe_fc,         'sothe_fc',           aoi, 250,          'Sothe et al. FC (kg/m2)');
  printStats(scanfi_fc,        'scanfi_fc',          aoi, EXPORT_SCALE, 'SCANFI v1.2 FC (kg/m2)');
  printStats(sothe_fc_unc,     'sothe_fc_unc',       aoi, 250,          'Sothe FC σ [per-pixel] (kg/m2)');
  printStats(forest_ens_mean,  'forest_ens_mean',    aoi, EXPORT_SCALE, 'Forest Equal-weight Mean [comparison] (kg/m2)');
  printStats(forest_wmean,     'forest_wmean',       aoi, EXPORT_SCALE, 'Forest Weighted Mean [primary] (kg/m2)');
  printStats(forest_wmean_sigma,'forest_wmean_sigma',aoi, EXPORT_SCALE, 'Forest Weighted σ (kg/m2)');
  print('  GEDI RF σ [constant]:  ' + GEDI_RF_SIGMA_KGPM2.toFixed(3) + ' kg/m²  (3-fold CV)');
  print('  SCANFI σ  [constant]:  ' + SCANFI_SIGMA_FC.toFixed(3) + ' kg/m²  (Guindon 2024)');
  print('');
  print('SOIL CARBON PRODUCTS');
  printStats(sothe_sc,         'sothe_sc',           aoi, 250, 'Sothe et al. SC (kg/m2)');
  printStats(sg_soc_1m,        'sg_soc_1m',          aoi, 250, 'SoilGrids OCS 0-100 cm (kg/m2)');
  printStats(sothe_sc_unc,     'sothe_sc_unc',       aoi, 250, 'Sothe SC σ [per-pixel] (kg/m2)');
  printStats(sg_sigma,         'sg_sigma',           aoi, 250, 'SoilGrids σ [derived] (kg/m2)');
  printStats(soil_ens_sd,      'soil_ens_sd',        aoi, 250, 'Inter-product Disagreement SD (kg/m2)');
  printStats(soil_ens_mean,    'soil_ens_mean',      aoi, 250, 'Soil Equal-weight Mean [comparison] (kg/m2)');
  printStats(soil_wmean,       'soil_wmean',         aoi, 250, 'Soil Weighted Mean [primary] (kg/m2)');
  printStats(soil_wmean_sigma, 'soil_wmean_sigma',   aoi, 250, 'Soil Weighted σ (kg/m2)');
  print('══════════════════════════════════════════════════════════════════');

  // ── Total carbon stocks ────────────────────────────────────────
  // AOI area in m² (geodesic, at export scale precision)
  var aoi_area_m2 = ee.Number(aoi.area(EXPORT_SCALE));
  var aoi_area_ha = aoi_area_m2.divide(1e4);

  // Reduce weighted means to get spatially averaged density
  var f_mean_ee = ee.Number(forest_wmean.updateMask(forest_nf_mask).reduceRegion({
    reducer: ee.Reducer.mean(), geometry: aoi,
    scale: EXPORT_SCALE, crs: EXPORT_CRS, maxPixels: 1e11, bestEffort: true
  }).get('forest_wmean'));

  var s_mean_ee = ee.Number(soil_wmean.reduceRegion({
    reducer: ee.Reducer.mean(), geometry: aoi,
    scale: 250, crs: EXPORT_CRS, maxPixels: 1e11, bestEffort: true
  }).get('soil_wmean'));

  // Total stock = mean density (kg/m²) × AOI area (m²) → kg → Mg (= t)
  var f_total_t  = f_mean_ee.multiply(aoi_area_m2).divide(1000);
  var s_total_t  = s_mean_ee.multiply(aoi_area_m2).divide(1000);
  var tot_total_t = f_total_t.add(s_total_t);

  ee.List([aoi_area_ha, f_mean_ee, f_total_t, s_mean_ee, s_total_t, tot_total_t])
    .evaluate(function(vals, err) {
      if (err || !vals) { print('Total C computation error: ' + err); return; }
      var ha     = vals[0];
      var f_den  = vals[1];   // kg/m²
      var f_tot  = vals[2];   // Mg = t C
      var s_den  = vals[3];
      var s_tot  = vals[4];
      var tot    = vals[5];
      var fmt    = function(v, d) { return v !== null ? v.toFixed(d) : 'null'; };
      var fmtInt = function(v) { return v !== null ? Math.round(v).toLocaleString() : 'null'; };

      print('');
      print('══ TOTAL CARBON STOCKS — WEIGHTED ENSEMBLE ══════════════════════');
      print('  AOI area:  ' + fmt(ha, 1) + ' ha');
      print('');
      print('  FOREST CARBON:');
      print('    Mean density:   ' + fmt(f_den*10, 1) + ' Mg C/ha  (' + fmt(f_den, 3) + ' kg/m²)');
      print('    Total stock:    ' + fmtInt(f_tot) + ' t C  (' + fmtInt(f_tot*1000) + ' kg C)');
      print('    [σ = ' + GEDI_RF_SIGMA_KGPM2.toFixed(3) + ' (GEDI RF)  |  per-pixel (Sothe)  |  ' +
            SCANFI_SIGMA_FC.toFixed(3) + ' (SCANFI)]');
      print('');
      print('  SOIL CARBON:');
      print('    Mean density:   ' + fmt(s_den*10, 1) + ' Mg C/ha  (' + fmt(s_den, 3) + ' kg/m²)');
      print('    Total stock:    ' + fmtInt(s_tot) + ' t C  (' + fmtInt(s_tot*1000) + ' kg C)');
      print('    [Sothe σ per-pixel  |  SoilGrids σ derived (standardized baseline + disagreement)]');
      print('');
      print('  TOTAL ECOSYSTEM:');
      print('    Mean density:   ' + fmt((f_den+s_den)*10, 1) + ' Mg C/ha  (' + fmt(f_den+s_den, 3) + ' kg/m²)');
      print('    Total stock:    ' + fmtInt(tot) + ' t C  (' + fmtInt(tot*1000) + ' kg C)');
      print('');
      print('  NOTE: Total stock = mean density × AOI area. This is a prior');
      print('  estimate from remote-sensing products. Field sampling will');
      print('  provide the calibrated stock with a defensible MOE and CI.');
      print('═══════════════════════════════════════════════════════════════════');
    });

  // ── Summary table (FeatureCollection → CSV export) ────────────
  var makeRow = function(source, desc, sigma_type, sigma_val, pool, res, st, scope) {
    return ee.Feature(null, {
      '1_source':      source,
      '2_description': desc,
      '2_carbon_scope': scope || 'AGB+BGB carbon',
      '3_sigma_type':  sigma_type,
      '3_sigma_mean':  sigma_val,
      '4_mean_kgm2':   st.mean,
      '4_sd_kgm2':     st.sd,
      '4_min_kgm2':    st.min,
      '4_max_kgm2':    st.max,
      'pool':          pool,
      'res_m':         res
    });
  };

  var gedi_rf_sigma_ee = ee.Number(GEDI_RF_SIGMA_KGPM2);
  var scanfi_sigma_ee  = ee.Number(SCANFI_SIGMA_FC);

  var tableRows = [
    // Forest
    makeRow('GEDI RF (this study)', 'Pilot→final 2-stage RF on GEDI L4A AGBD; AGB×' + AGB_TO_C_KGM2.toFixed(3) + '; 3-fold CV σ',
      '3-fold CV RMSE × ' + AGB_TO_C_KGM2.toFixed(3), gedi_rf_sigma_ee, 'Forest Carbon', EXPORT_SCALE,
      getStats(forest_rf_pred, 'forest_carbon_kgm2', EXPORT_SCALE), 'AGB+BGB carbon (from AGBD)'),
    makeRow('Sothe et al. FC', 'McMaster/WWF-Canada national FC map (kg/m2)',
      'Per-pixel product unc', getStats(sothe_fc_unc,'sothe_fc_unc',250).mean, 'Forest Carbon', 250,
      getStats(sothe_fc, 'sothe_fc', 250), 'AGB+BGB+dead carbon (native)'),
    makeRow('SCANFI v1.2 FC', 'AGB × ' + AGB_TO_C_KGM2.toFixed(3) + ' → kg C/m². Guindon 2024 RMSE 38.70 t/ha',
      'Constant (published RMSE)', scanfi_sigma_ee, 'Forest Carbon', 25,
      getStats(scanfi_fc, 'scanfi_fc', EXPORT_SCALE), 'AGB+BGB carbon (from AGB)'),
    makeRow('Forest Weighted Mean (this study)',
      'Inv-var weighted: GEDI RF + Sothe FC + SCANFI. w_i = 1/σ_i²',
      'Inv-var combined σ', getStats(forest_wmean_sigma,'forest_wmean_sigma',EXPORT_SCALE).mean,
      'Forest Carbon', EXPORT_SCALE, getStats(forest_wmean, 'forest_wmean', EXPORT_SCALE), 'AGB+BGB carbon'),
    makeRow('Forest Equal-weight Mean [comparison]', 'GEDI RF + Sothe FC + SCANFI, equal weights',
      'Ensemble SD', getStats(forest_ens_sd,'forest_ens_sd',EXPORT_SCALE).mean,
      'Forest Carbon [comparison]', EXPORT_SCALE, getStats(forest_ens_mean,'forest_ens_mean',EXPORT_SCALE), 'AGB+BGB carbon'),
    // Soil
    makeRow('Sothe et al. SC', 'McMaster/WWF-Canada national SC map (kg/m2)',
      'Per-pixel product unc', getStats(sothe_sc_unc,'sothe_sc_unc',250).mean, 'Soil Carbon', 250,
      getStats(sothe_sc, 'sothe_sc', 250), 'Soil organic carbon 0-1 m'),
    makeRow('SoilGrids OCS 0-100 cm', 'Depth-integrated OCS from SoilGrids v2.0 (kg/m2)',
      'Derived (std Sothe unc + inter-product SD)',
      getStats(sg_sigma,'sg_sigma',250).mean, 'Soil Carbon', 250,
      getStats(sg_soc_1m, 'sg_soc_1m', 250), 'Soil organic carbon 0-1 m'),
    makeRow('Soil Weighted Mean (this study)',
      'Inv-var weighted: Sothe SC + SoilGrids. SoilGrids σ derived.',
      'Inv-var combined σ', getStats(soil_wmean_sigma,'soil_wmean_sigma',250).mean,
      'Soil Carbon', 250, getStats(soil_wmean, 'soil_wmean', 250), 'Soil organic carbon 0-1 m'),
    makeRow('Soil Equal-weight Mean [comparison]', 'Sothe SC + SoilGrids, equal weights',
      'Ensemble SD (inter-product disagree)', getStats(soil_ens_sd,'soil_ens_sd',250).mean,
      'Soil Carbon [comparison]', 250, getStats(soil_ens_mean,'soil_ens_mean',250), 'Soil organic carbon 0-1 m'),
    // Total
    makeRow('Total Ecosystem Carbon — Weighted', 'Forest Weighted Mean + Soil Weighted Mean',
      'Propagated inv-var', null, 'Total Ecosystem', EXPORT_SCALE,
      getStats(total_ecosystem_c_w, 'total_ec_weighted', EXPORT_SCALE), 'Forest + soil carbon')
  ];

  var summaryTable = ee.FeatureCollection(tableRows.filter(function(r){return r!==null;}));

  // ── Exports ────────────────────────────────────────────────────
  var exportRaster = function(img, desc) {
    Export.image.toDrive({ image: img.toFloat(), description: desc,
      folder: EXPORT_FOLDER, region: aoi, scale: EXPORT_SCALE, maxPixels: 1e13, crs: EXPORT_CRS });
  };

  // Weighted ensemble outputs (primary)
  exportRaster(forest_wmean,        'Forest_Weighted_Mean_kgm2');
  exportRaster(forest_wmean_sigma,  'Forest_Weighted_Sigma_kgm2');
  exportRaster(soil_wmean,          'Soil_Weighted_Mean_kgm2');
  exportRaster(soil_wmean_sigma,    'Soil_Weighted_Sigma_kgm2');
  exportRaster(sg_sigma,            'SoilGrids_Derived_Sigma_kgm2');
  exportRaster(total_ecosystem_c_w, 'Total_Ecosystem_Carbon_Weighted_kgm2');

  // Equal-weight ensemble (comparison)
  exportRaster(forest_ens_mean,     'Forest_EqualWeight_Mean_kgm2');
  exportRaster(forest_ens_sd,       'Forest_EqualWeight_SD_kgm2');
  exportRaster(soil_ens_mean,       'Soil_EqualWeight_Mean_kgm2');
  exportRaster(soil_ens_sd,         'Soil_InterProduct_Disagree_SD_kgm2');
  exportRaster(total_ecosystem_c,   'Total_Ecosystem_Carbon_EqualWeight_kgm2');

  // Uncertainty
  exportRaster(sothe_fc_unc,        'Sothe_FC_Uncertainty_kgm2');
  exportRaster(sothe_sc_unc,        'Sothe_SC_Uncertainty_kgm2');
  exportRaster(forest_uncertainty,  'Forest_RSS_Uncertainty_kgm2');
  exportRaster(soil_uncertainty,    'Soil_RSS_Uncertainty_kgm2');

  // Source products
  exportRaster(forest_rf_pred,      'GEDI_RF_Forest_Carbon_kgm2');
  exportRaster(sothe_fc,            'Sothe_FC_kgm2');
  exportRaster(scanfi_fc,           'SCANFI_Forest_Carbon_kgm2');
  exportRaster(sothe_sc,            'Sothe_SC_kgm2');
  exportRaster(sg_soc_1m,           'SoilGrids_OCS_0_100cm_kgm2');

  // Data source maps
  exportRaster(dsm_forest,          'Forest_Data_Source_Map');
  exportRaster(dsm_soil,            'Soil_Data_Source_Map');

  // Tables
  Export.table.toDrive({ collection: summaryTable, description: 'Carbon_Summary_Table_v4_6',
    folder: EXPORT_FOLDER, fileFormat: 'CSV' });

  if (forest_sampling_pts) {
    Export.table.toDrive({ collection: forest_sampling_pts,
      description: 'Forest_Sampling_Points_Neyman_KML', folder: EXPORT_FOLDER, fileFormat: 'KML' });
    Export.table.toDrive({ collection: forest_sampling_pts,
      description: 'Forest_Sampling_Points_Neyman_CSV', folder: EXPORT_FOLDER, fileFormat: 'CSV' });
  } else print('Forest sampling points not generated - run Step 8.');

  if (soil_sampling_pts) {
    Export.table.toDrive({ collection: soil_sampling_pts,
      description: 'Soil_Sampling_Points_Neyman_KML', folder: EXPORT_FOLDER, fileFormat: 'KML' });
    Export.table.toDrive({ collection: soil_sampling_pts,
      description: 'Soil_Sampling_Points_Neyman_CSV', folder: EXPORT_FOLDER, fileFormat: 'CSV' });
  } else print('Soil sampling points not generated - run Step 8.');

  markDone(8);
  setStatus('Step 9 complete - check Tasks panel to run all exports.');
  if (typeof onDone === 'function') onDone();
}


// ─────────────────────────────────────────────────────────────────
// SECTION G — PARTNER ORCHESTRATION (asset preflight + Run All)
// ─────────────────────────────────────────────────────────────────

// [0] Verify that every asset the pipeline reads is accessible to this
// account BEFORE running anything, so a partner gets a clear, early
// diagnosis instead of a mid-run permission error.
function checkAssetAccess() {
  setStatus('[0] Checking Earth Engine asset access...');
  print('');
  print('══ ASSET ACCESS PREFLIGHT ════════════════════════════════════════');
  print('  FAIL on a REQUIRED asset → ask the owner to share it with your');
  print('  account, or replace the ID in SECTION A with one you can read.');
  print('  (Public datasets — GEDI, SoilGrids, SCANFI, sat-io — are not');
  print('  re-checked here; only project-private assets are listed.)');
  print('');

  var checks = [
    { id: AOI_ASSET,          type: 'table', required: true,  label: 'AOI boundary' },
    { id: SOTHE_FC_UNC_ASSET, type: 'image', required: true,  label: 'Sothe FC uncertainty' },
    { id: SOTHE_SC_UNC_ASSET, type: 'image', required: true,  label: 'Sothe SC uncertainty' },
    { id: WOSIS_ASSET,        type: 'table', required: false, label: 'WOSIS field data (Step 1 only)' },
    { id: CANPEAT_ASSET,      type: 'table', required: false, label: 'CanPeat field data (Step 1 only)' },
    { id: COMBINED_ASSET,     type: 'table', required: false, label: 'Combined profiles (Step 4b only)' }
  ];

  var remaining = checks.length, reqFail = 0;
  var finalize = function() {
    remaining -= 1;
    if (remaining > 0) return;
    print('');
    if (reqFail === 0) {
      print('  ✓ All REQUIRED assets readable — safe to run the pipeline.');
      setStatus('[0] Access OK — click [▶ RUN ALL] or Step [2].');
    } else {
      print('  ✗ ' + reqFail + ' REQUIRED asset(s) NOT readable — fix sharing first.');
      setStatus('[0] ' + reqFail + ' required asset(s) inaccessible — see console.');
    }
    print('══════════════════════════════════════════════════════════════════');
  };

  checks.forEach(function(c) {
    var probe = (c.type === 'image')
      ? ee.Image(c.id).bandNames()
      : ee.FeatureCollection(c.id).limit(0).size();
    probe.evaluate(function(res, err) {
      var ok = !err && res !== null && res !== undefined;
      print('  [' + (ok ? '✓ PASS' : '✗ FAIL') + '] ' +
            (c.required ? 'REQUIRED ' : 'optional ') + c.label);
      print('        ' + c.id + (err ? '\n        → ' + err : ''));
      if (!ok && c.required) reqFail += 1;
      finalize();
    });
  });
}

// [▶ RUN ALL] Chain Steps 2→9. Each step calls the next from its own
// completion callback (set up via the onDone parameter), so the async
// GEDI training / 3-fold CV / soil ensemble all finish before the next
// step starts. Step 1 (field data) and [4b] snapshot are optional and
// are intentionally NOT run here.
function runAll() {
  setStatus('▶ RUN ALL: starting Steps 2 → 9...');
  print('');
  print('▶ RUN ALL — executing Steps 2 through 9 in sequence. This can take');
  print('  several minutes; watch the PROGRESS panel. When it finishes, open');
  print('  the Tasks tab (top-right) to run the exports.');
  step2_importRasterPriors(function() {
  step3_importGEDI(function() {
  step4_buildCovariates(function() {
  step5_trainFinalGEDI(function() {
  step6_buildSoilModel(function() {
  step7_buildEnsemble(function() {
  step8_generateSampling(function() {
  step9_reportsAndExports(function() {
    setStatus('▶ RUN ALL complete — open the Tasks tab to run exports.');
    print('▶ RUN ALL complete. Open the Tasks tab (top-right) to run exports.');
  }); }); }); }); }); }); }); });
}


// ─────────────────────────────────────────────────────────────────
// BOOT
// ─────────────────────────────────────────────────────────────────
print("=== Charlie's Place KBA - Forest Carbon Assessment v4.6 ===");
print('CRS: ' + EXPORT_CRS + ' | Scale: ' + EXPORT_SCALE + ' m');
print('Forest: ' + N_FOREST_SAMPLES + ' pts | Soil: ' + N_SOIL_SAMPLES + ' pts (EXACT)');
print('Carbon: AGB→C factor = ' + AGB_TO_C_KGM2.toFixed(3) +
      '  (×(1+' + BGB_RATIO + ') BGB ×' + CARBON_FRACTION + ' C-frac ×' + MGHA_TO_KGM2 + ' unit)');
print('Ensemble: inverse-variance weighting (w_i = 1/σ_i²) per pool');
print('Forest σ: GEDI RF = 3-fold CV | Sothe FC = per-pixel | SCANFI = ' + SCANFI_SIGMA_FC.toFixed(3) + ' kg/m²');
print('Soil   σ: Sothe SC = per-pixel | SoilGrids = standardized(Sothe unc + inter-product SD)');
print('Sampling: Neyman (largest-remainder, exact n) | LC: ESA WorldCover v200 × ' + N_UNC_BINS + ' unc bins');
print('Partner flow: [0] Check Asset Access → [▶ RUN ALL] → Tasks tab → run exports.');

try { initUI(); } catch(e) {
  print('Panel already rendered - hard refresh (Ctrl+Shift+R) to fully reset.');
}
