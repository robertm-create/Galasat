// ============================================================
// GalaSat v1.0 — trainingPoints architecture, 506-line build
//
// MEMORY TAG: GalaSat_v1.0_trainingPoints_506
// SOURCE: claude.ai chat "Seasonal flood zone mapping in Earth Engine"
//         https://claude.ai/chat/706da8c8-d698-4916-9a82-9da882bced38
//         pasted attachment, 506 lines
// EXTRACTED: 2026-04-25
//
// LINEAGE NOTE: Parallel to the SuperSat lineage. Trains on the 340-point
// asset directly (not the polygon assets used by V0–V3). 86.7% overall
// accuracy per legend. Self-labeled "GalaSat v1.0".
//
// Architecture:
//   - AOI = trainingPoints.geometry().bounds().buffer(20000) — dynamic, not hardcoded
//   - 16-band feature stack (B2–B12 + NDVI/NDWI/MNDWI/NDTI/BSI/IronOxide)
//   - Binary RF (mining vs clean vegetation, NDVI > 0.6 reference points)
//   - Dual composite: dry season disturbance + wet season turbid water
//   - 4 independent water raster layers (all waterways, tier1, adjacent, exposed communities)
//   - 12-year time series (2014–2025) using L8/L9 merged
//   - F1 score per class (first appearance in lineage; AOI_WATER_LAYER also has it)
//
// NOTE: Whitespace partially collapsed during DOM extraction. Logic and
// semicolons intact — runnable.

// ── ASSETS ─────────────────────────────────────────────────────────────
var trainingPoints = ee.FeatureCollection(
  'projects/galamsey-monotoring/assets/infn8vzn_mercury_training_340'
);
var AOI = trainingPoints.geometry().bounds().buffer(20000);

// ── OSM ASSETS ─────────────────────────────────────────────────────────
var roads = ee.FeatureCollection(
  'projects/galamsey-monotoring/assets/ghana_osm_roads'
).filterBounds(AOI)
 .filter(ee.Filter.inList('fclass', [
   'motorway','trunk','primary','secondary'
 ]));
var waterways = ee.FeatureCollection(
  'projects/galamsey-monotoring/assets/ghana_osm_waterways'
).filterBounds(AOI);
var places = ee.FeatureCollection(
  'projects/galamsey-monotoring/assets/ghana_osm_places'
).filterBounds(AOI);

// ── CONFIRMED CONTAMINATED RIVERS ──────────────────────────────────────
var tier1Rivers = waterways.filter(
  ee.Filter.inList('name', [
    'Pra','Pra River',
    'Offin','Offin River',
    'Ankobra','Ankobra River',
    'Birim','Birim River',
    'Tano','Tano River',
    'Densu','Densu River',
    'Oda','Oda River',
    'Bonsa','Bonsa River',
    'Oti','Oti River',
    'Fena','Fena River',
    'Subin','Subin River',
    'Anum','Anum River',
    'Afram','Afram River',
    'Anikoko','Bodwire',
    'Asesree','Assaman'
  ])
);

// ── GLOBAL MINING WATCH ────────────────────────────────────────────────
var globalMines = ee.FeatureCollection(
  'projects/sat-io/open-datasets/global-mining/global_mining_polygons'
).filterBounds(AOI);

var minesClipped = globalMines.map(function(f) {
  return f.intersection(AOI, ee.ErrorMargin(1));
});

var totalMinesHa = ee.Number(minesClipped.map(function(f) {
  return f.set('area_ha', f.geometry().area(1).divide(10000));
}).aggregate_sum('area_ha'));

var miningBuffer = globalMines.map(function(f) {
  return f.buffer(2000);
});
var miningBufferImage = ee.Image(0).byte().paint(miningBuffer, 1);

print('── Mining Footprints ──────────────────');
print('Mining polygons in AOI:', globalMines.size());
print('Mining footprint area (hectares):', totalMinesHa);
print('Confirmed river segments in AOI:', tier1Rivers.size());
print('Total named places in AOI:', places.size());

// ── CLOUD MASK S2 ──────────────────────────────────────────────────────
function maskS2clouds(image) {
  var qa = image.select('QA60');
  var cloudBitMask = 1 << 10;
  var cirrusBitMask = 1 << 11;
  var mask = qa.bitwiseAnd(cloudBitMask).eq(0)
    .and(qa.bitwiseAnd(cirrusBitMask).eq(0));
  return image.updateMask(mask).divide(10000);
}

// ── DRY SEASON 2025-26 — disturbance detection ─────────────────────────
var s2dry = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
  .filterBounds(AOI)
  .filterDate('2025-11-01', '2026-02-28')
  .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 20))
  .map(maskS2clouds)
  .median()
  .clip(AOI);

// ── WET SEASON 2025 — turbid water detection ───────────────────────────
var s2wet2025 = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
  .filterBounds(AOI)
  .filterDate('2025-06-01', '2025-09-30')
  .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 30))
  .map(maskS2clouds)
  .median()
  .clip(AOI);

print('── Wet Season 2025 image count ────────');
print('Images in 2025 wet season composite:',
  ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
    .filterBounds(AOI)
    .filterDate('2025-06-01', '2025-09-30')
    .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 30))
    .size()
);

// ── DRY SEASON INDICES ─────────────────────────────────────────────────
var NDVI  = s2dry.normalizedDifference(['B8','B4']).rename('NDVI');
var NDWI  = s2dry.normalizedDifference(['B3','B8']).rename('NDWI');
var MNDWI = s2dry.normalizedDifference(['B3','B11']).rename('MNDWI');
var NDTI  = s2dry.normalizedDifference(['B4','B3']).rename('NDTI');
var BSI   = s2dry.expression(
  '((SWIR + RED) - (NIR + BLUE)) / ((SWIR + RED) + (NIR + BLUE))',
  {SWIR:s2dry.select('B11'), RED:s2dry.select('B4'),
   NIR:s2dry.select('B8'), BLUE:s2dry.select('B2')}
).rename('BSI');
var IOR = s2dry.select('B4').divide(s2dry.select('B2')).rename('IronOxide');

// ── WET SEASON 2025 INDICES ────────────────────────────────────────────
var MNDWI_wet = s2wet2025.normalizedDifference(['B3','B11']).rename('MNDWI_wet');
var NDTI_wet  = s2wet2025.normalizedDifference(['B4','B3']).rename('NDTI_wet');
var IOR_wet   = s2wet2025.select('B4').divide(s2wet2025.select('B2')).rename('IOR_wet');
var NDVI_wet  = s2wet2025.normalizedDifference(['B8','B4']).rename('NDVI_wet');

// ── FEATURE STACK ──────────────────────────────────────────────────────
var imageStack = s2dry.select(['B2','B3','B4','B5','B6','B7','B8','B8A','B11','B12'])
  .addBands(NDVI).addBands(NDWI).addBands(MNDWI)
  .addBands(NDTI).addBands(BSI).addBands(IOR);

// ── CLEAN REFERENCE POINTS ─────────────────────────────────────────────
var cleanForestMask = NDVI.gt(0.6);
var cleanPoints = imageStack.updateMask(cleanForestMask).sample({
  region: AOI, scale: 30, numPixels: 278, seed: 42, geometries: true
});
cleanPoints = cleanPoints.map(function(f) {
  return f.set('class', 0);
});

// ── SAMPLE MERCURY POINTS ──────────────────────────────────────────────
var mercurySampled = imageStack.sampleRegions({
  collection: trainingPoints.filter(ee.Filter.eq('class', 1)),
  properties: ['class'], scale: 10, tileScale: 8
});

// ── MERGE AND SPLIT ────────────────────────────────────────────────────
var allSamples = mercurySampled.merge(cleanPoints)
  .randomColumn('random', 42);
var trainSet = allSamples.filter(ee.Filter.lt('random', 0.8));
var valSet   = allSamples.filter(ee.Filter.gte('random', 0.8));

print('── STEP 7: ACCURACY CHECK ─────────────');
print('Total samples:', allSamples.size());
print('Training set:', trainSet.size());
print('Validation set:', valSet.size());

// ── TRAIN RANDOM FOREST ────────────────────────────────────────────────
var bands = ['B2','B3','B4','B5','B6','B7','B8','B8A','B11','B12',
             'NDVI','NDWI','MNDWI','NDTI','BSI','IronOxide'];

var classifier = ee.Classifier.smileRandomForest({
  numberOfTrees: 200, seed: 42
}).train({
  features: trainSet, classProperty: 'class', inputProperties: bands
});

// ── CLASSIFY ───────────────────────────────────────────────────────────
var classified = imageStack.classify(classifier);

// ── MASKS ──────────────────────────────────────────────────────────────
var worldCover = ee.ImageCollection('ESA/WorldCover/v200').first().clip(AOI);
var urbanMask = worldCover.neq(50).and(worldCover.neq(40));
var roadBuffered = roads.map(function(f) { return f.buffer(50); });
var roadMask = ee.Image(1).byte().paint(roadBuffered, 0);
var finalMask = urbanMask.and(roadMask);

// ── CONTAMINATION ONLY ─────────────────────────────────────────────────
var contaminationOnly = classified
  .updateMask(finalMask)
  .updateMask(classified.eq(1))
  .selfMask();

// ── FOREST LAYER ───────────────────────────────────────────────────────
var forest = worldCover.eq(10).selfMask();

// ── TURBID WATER — WET SEASON 2025 — FIXED SCALE 30 ─────────────────────
var turbidWater2025 = MNDWI_wet.gt(0)
  .and(NDTI_wet.gt(0.05))
  .and(IOR_wet.gt(1.05))
  .and(NDVI_wet.lt(0.3))
  .selfMask()
  .rename('TurbidWater');

// ── FOUR INDEPENDENT RASTER WATER LAYERS ───────────────────────────────
// LAYER 1 — BLUE — All waterways — width 1
var allWaterwaysRaster = ee.Image().byte().paint({
  featureCollection: waterways, color: 1, width: 1
});

// LAYER 2 — CRIMSON — Field-verified named rivers — width 3
var tier1Raster = ee.Image().byte().paint({
  featureCollection: tier1Rivers, color: 1, width: 3
});

// LAYER 3 — PURPLE — Waterways adjacent to confirmed mining
var adjacentWaterwaysRaster = ee.Image().byte().paint(
  waterways.map(function(f) { return f.buffer(30); }), 1
).updateMask(miningBufferImage);

// LAYER 4 — MAGENTA — Exposed communities — 100m dots
var placesWithExposure = places.map(function(f) {
  var buf = f.buffer(5000);
  var disturbanceInBuffer = contaminationOnly.reduceRegion({
    reducer: ee.Reducer.sum(),
    geometry: buf.geometry(),
    scale: 100, maxPixels: 1e9
  });
  var hasDisturbance = ee.Number(
    disturbanceInBuffer.get('classification')).gt(0);
  return f.set('exposed', hasDisturbance);
});
var exposedCommunities = placesWithExposure
  .filter(ee.Filter.eq('exposed', 1));
var exposedCommunitiesRaster = ee.Image().byte().paint({
  featureCollection: exposedCommunities.map(function(f) {
    return f.buffer(100);
  }),
  color: 1
});

// ── AREA CALCULATIONS — FIXED SCALE 30 ─────────────────────────────────
var turbidAreaDict = turbidWater2025
  .multiply(ee.Image.pixelArea()).divide(10000)
  .reduceRegion({
    reducer: ee.Reducer.sum(),
    geometry: AOI, scale: 30, maxPixels: 1e13
  });
print('── STEP 2: Turbid Water (wet season 2025) ─');
print('Turbid Water Area (hectares):', turbidAreaDict.get('TurbidWater'));

var contaminationAreaDict = contaminationOnly
  .multiply(ee.Image.pixelArea()).divide(10000)
  .reduceRegion({
    reducer: ee.Reducer.sum(),
    geometry: AOI, scale: 30, maxPixels: 1e13
  });
print('── STEP 3: Mining Disturbance 2026 ────');
print('Mining Disturbance Area (hectares):', contaminationAreaDict.get('classification'));

print('── COMMUNITY EXPOSURE ─────────────────');
print('Total communities in AOI:', places.size());
print('Communities within 5km of disturbance:', exposedCommunities.size());

// ── VALIDATION ─────────────────────────────────────────────────────────
var validated = valSet.classify(classifier);
var errorMatrix = validated.errorMatrix('class', 'classification');
print('── Confusion Matrix ───────────────────');
print('Confusion Matrix:', errorMatrix);
print('Overall Accuracy:', errorMatrix.accuracy());
print('Kappa:', errorMatrix.kappa());
print('Producers Accuracy:', errorMatrix.producersAccuracy());
print('Consumers Accuracy:', errorMatrix.consumersAccuracy());
print('F1 per class:', errorMatrix.fscore());

// ── STEP 6: TIME SERIES 2014-2025 ──────────────────────────────────────
function getYearTurbidHa(year) {
  var startDate = year + '-06-01';
  var endDate = year + '-09-30';
  function maskL8(image) {
    var qa = image.select('QA_PIXEL');
    var mask = qa.bitwiseAnd(1 << 3).eq(0)
      .and(qa.bitwiseAnd(1 << 4).eq(0));
    return image.updateMask(mask)
      .select(['SR_B3','SR_B4','SR_B6'])
      .multiply(0.0000275).add(-0.2);
  }
  function maskL9(image) {
    var qa = image.select('QA_PIXEL');
    var mask = qa.bitwiseAnd(1 << 3).eq(0)
      .and(qa.bitwiseAnd(1 << 4).eq(0));
    return image.updateMask(mask)
      .select(['SR_B3','SR_B4','SR_B6'])
      .multiply(0.0000275).add(-0.2);
  }
  var L8 = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
    .filterBounds(AOI).filterDate(startDate, endDate)
    .filter(ee.Filter.lt('CLOUD_COVER', 30)).map(maskL8);
  var L9 = ee.ImageCollection('LANDSAT/LC09/C02/T1_L2')
    .filterBounds(AOI).filterDate(startDate, endDate)
    .filter(ee.Filter.lt('CLOUD_COVER', 30)).map(maskL9);
  var merged = L8.merge(L9);
  var result = ee.Algorithms.If(
    merged.size().gte(3),
    ee.Number(
      merged.median().clip(AOI)
        .normalizedDifference(['SR_B3','SR_B6']).rename('MNDWI')
        .gt(0)
        .and(
          merged.median().clip(AOI)
            .normalizedDifference(['SR_B4','SR_B3']).rename('NDTI')
            .gt(0.05)
        )
        .rename('turbid_ha')
        .multiply(ee.Image.pixelArea()).divide(10000)
        .reduceRegion({
          reducer: ee.Reducer.sum(), geometry: AOI,
          scale: 30, maxPixels: 1e13, bestEffort: true
        }).get('turbid_ha')
    ),
    null
  );
  return result;
}

var years = [2014,2015,2016,2017,2018,2019,2020,2021,2022,2023,2024,2025];
var featureList = years.map(function(year) {
  var ha = getYearTurbidHa(year);
  return ee.Feature(null, {year: year, turbid_ha: ha});
});
var timeSeriesFC = ee.FeatureCollection(featureList);

var chart = ui.Chart.feature.byFeature({
  features: timeSeriesFC, xProperty: 'year', yProperties: ['turbid_ha']
})
.setChartType('LineChart')
.setOptions({
  title: 'GalaSat: Mining-Affected Turbid Water 2014-2025 — Wet Season (Ashanti Pilot)',
  hAxis: {title: 'Year', format: '####'},
  vAxis: {title: 'Turbid Water Area (hectares)'},
  colors: ['#8B0000'], lineWidth: 3, pointSize: 6,
  backgroundColor: 'white', legend: {position: 'none'},
  interpolateNulls: false
});
print('── STEP 6: Time Series Chart ───────────');
print(chart);

// ── VISUALIZE ─────────────────────────────────────────────────────────
Map.centerObject(AOI, 10);
var aoiOutline = ee.Image().byte().paint({
  featureCollection: ee.FeatureCollection([ee.Feature(AOI)]),
  color: 1, width: 2
});

Map.addLayer(s2dry, {bands: ['B4','B3','B2'], min: 0, max: 0.3}, 'True Colour 2025-26');
Map.addLayer(forest, {palette: ['1A5C1A']}, 'Forest Cover');
Map.addLayer(contaminationOnly, {palette: ['FF0000'], opacity: 0.6}, 'Mining Disturbance Signal');
Map.addLayer(turbidWater2025, {palette: ['FF8C00'], opacity: 0.8}, 'Turbid Water (wet season 2025)');
Map.addLayer(globalMines, {color: 'FFFF00'}, 'Mining Footprints — Global Mining Watch');
Map.addLayer(roads, {color: '888888'}, 'Major Roads');
Map.addLayer(allWaterwaysRaster, {palette: ['4FC3F7'], opacity: 0.5}, 'All Waterways (blue)');
Map.addLayer(tier1Raster, {palette: ['8B0000']}, 'Field-Verified Contaminated Rivers (crimson)');
Map.addLayer(adjacentWaterwaysRaster, {palette: ['9C27B0']}, 'Waterways Adjacent to Confirmed Mining (purple)');
Map.addLayer(exposedCommunitiesRaster, {palette: ['FF00FF']}, 'Communities within 5km (magenta)');
Map.addLayer(trainingPoints, {color: 'FF6600'}, 'Training Points');
Map.addLayer(aoiOutline, {palette: ['FFFFFF'], opacity: 0.8}, 'Pilot AOI Boundary');

// ── LEGEND ─────────────────────────────────────────────────────────────
var legend = ui.Panel({
  style: {
    position: 'bottom-left', padding: '8px 15px', backgroundColor: 'white'
  }
});
legend.add(ui.Label({
  value: 'GalaSat v1.0 — INFN8 VZN',
  style: {fontWeight: 'bold', fontSize: '14px', margin: '0 0 2px 0', color: '#8B0000'}
}));
legend.add(ui.Label({
  value: 'Ashanti Pilot, Ghana | March 2026',
  style: {fontSize: '11px', margin: '0 0 6px 0', color: '#555555'}
}));
var makeRow = function(color, label) {
  var colorBox = ui.Label({
    style: {backgroundColor: color, padding: '8px', margin: '0 6px 4px 0'}
  });
  var description = ui.Label({
    value: label,
    style: {margin: '0 0 4px 6px', fontSize: '12px'}
  });
  return ui.Panel({
    widgets: [colorBox, description],
    layout: ui.Panel.Layout.Flow('horizontal')
  });
};
legend.add(makeRow('#1A5C1A', 'Forest Cover'));
legend.add(makeRow('#FF0000', 'Mining Disturbance (Nov 2025-Feb 2026)'));
legend.add(makeRow('#FF8C00', 'Turbid Water (Jun-Sep 2025)'));
legend.add(makeRow('#FFFF00', 'Mining Footprints (Global Mining Watch)'));
legend.add(makeRow('#888888', 'Major Roads'));
legend.add(makeRow('#4FC3F7', 'All Waterways'));
legend.add(makeRow('#8B0000', 'Field-Verified Contaminated Rivers'));
legend.add(makeRow('#9C27B0', 'Waterways Adjacent to Confirmed Mining'));
legend.add(makeRow('#FF00FF', 'Communities within 5km'));
legend.add(makeRow('#FF6600', 'Training Points (340)'));
legend.add(makeRow('#FFFFFF', 'Pilot AOI Boundary'));
legend.add(ui.Label({
  value: 'Sentinel-2 SR | Dual composite 2025-26 | 86.7% OA | EPSG:4326',
  style: {fontSize: '10px', color: '#555555', margin: '6px 0 0 0'}
}));
legend.add(ui.Label({
  value: 'Sources: Joy News | Al Jazeera | PMC 2024 | GWCL | Global Mining Watch',
  style: {fontSize: '10px', color: '#555555', margin: '2px 0 0 0'}
}));
Map.add(legend);

// ── EXPORTS ────────────────────────────────────────────────────────────
Export.image.toDrive({
  image: contaminationOnly,
  description: 'INFN8VZN_Mercury_Disturbance_2026',
  folder: 'GalaSat',
  fileNamePrefix: 'infn8vzn_mercury_disturbance_2026',
  region: AOI, scale: 10, crs: 'EPSG:4326', maxPixels: 1e13, fileFormat: 'GeoTIFF'
});
Export.image.toDrive({
  image: turbidWater2025,
  description: 'INFN8VZN_Turbid_Water_2025',
  folder: 'GalaSat',
  fileNamePrefix: 'infn8vzn_turbid_water_2025',
  region: AOI, scale: 10, crs: 'EPSG:4326', maxPixels: 1e13, fileFormat: 'GeoTIFF'
});
print('── Exports submitted — check Tasks tab ──');
