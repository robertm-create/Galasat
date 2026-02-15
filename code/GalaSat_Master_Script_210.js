// =============================================================================
// GALASAT MASTER SCRIPT — INFN8 VZN
// Layers: Contamination Classification | Turbid Water | Flood Safety
// Basins: Pra, Ankobra, Birim, Tano, Offin | AOI: Ashanti Region, Ghana
// =============================================================================
//
// MEMORY TAG: GalaSat_Master_Script_210
// SOURCE: claude.ai chat "Seasonal flood zone mapping in Earth Engine"
//         https://claude.ai/chat/706da8c8-d698-4916-9a82-9da882bced38
//         pasted attachment, 210 lines
// EXTRACTED: 2026-04-25
//
// LINEAGE NOTE: Different architecture from the trainingPoints branch siblings
// in this same chat. Hardcoded AOI rectangle (not dynamic from polygon bounds),
// single 2024-2025 composite, RF set to PROBABILITY mode with > 0.75 threshold,
// balanced class 1/class 0 sampling (limit class 0 to match class 1 size),
// dry/wet JRC season layers shown separately. Likely an earlier "master script"
// draft from before the trainingPoints branch + later builds diverged.
//
// CODE QUIRK: There are two declarations of `var contaminationOnly` and `var worldCover`
// in this file. The first definition uses probability > 0.75 + finalMask; the
// second redeclares using the simpler `classified.updateMask(finalMask)` pattern.
// The second declaration shadows the first. Captured as-pasted; clean up if used.
//
// NOTE: Whitespace partially collapsed during DOM extraction.

Map.setCenter(-1.8, 6.3, 9);
Map.setOptions('SATELLITE');

// ── AOI ────────────────────────────────────────────────────────────────
var AOI = ee.Geometry.Rectangle([-3.00, 5.80, -0.90, 7.60], null, false);

// ── SENTINEL-2 COMPOSITE ───────────────────────────────────────────────
function maskS2scl(img) {
  var scl = img.select('SCL');
  var mask = scl.neq(3)   // cloud shadow
    .and(scl.neq(8))      // cloud medium probability
    .and(scl.neq(9))      // cloud high probability
    .and(scl.neq(10));    // thin cirrus
  return img.updateMask(mask);
}

var s2 = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
  .filterBounds(AOI)
  .filterDate('2024-01-01', '2025-01-01')
  .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 20))
  .map(maskS2scl)
  .median()
  .clip(AOI);

// ── SPECTRAL INDICES ───────────────────────────────────────────────────
var NDVI  = s2.normalizedDifference(['B8','B4']).rename('NDVI');
var NDWI  = s2.normalizedDifference(['B3','B8']).rename('NDWI');
var MNDWI = s2.normalizedDifference(['B3','B11']).rename('MNDWI');
var NDTI  = s2.normalizedDifference(['B4','B3']).rename('NDTI');
var BSI   = s2.expression(
  '((SWIR + RED) - (NIR + BLUE)) / ((SWIR + RED) + (NIR + BLUE))',
  {SWIR:s2.select('B11'), RED:s2.select('B4'),
   NIR:s2.select('B8'), BLUE:s2.select('B2')}
).rename('BSI');
var IOR = s2.select('B4').divide(s2.select('B2')).rename('IronOxide');

var imageStack = s2.select(['B2','B3','B4','B5','B6','B7','B8','B8A','B11','B12'])
  .addBands(NDVI).addBands(NDWI).addBands(MNDWI)
  .addBands(NDTI).addBands(BSI).addBands(IOR);

// ── STEP 2: TURBID WATER (water pixels only) ───────────────────────────
// Water mask: MNDWI > 0 = confirmed water surface
// Turbidity: NDTI > 0 on water pixels only = mercury-contaminated water
var waterMask = MNDWI.gt(0);
var turbidWater = MNDWI.gt(0).and(NDTI.gt(0))
  .selfMask()
  .rename('TurbidWater');

var turbidAreaDict = ee.Image.pixelArea()
  .updateMask(turbidWater)
  .rename('area')
  .divide(10000)
  .reduceRegion({
    reducer: ee.Reducer.sum(),
    geometry: AOI,
    scale: 20, maxPixels: 1e13, bestEffort: true
  });

// ── STEP 3: CONTAMINATION CLASSIFIER ───────────────────────────────────
var trainingPoints = ee.FeatureCollection('projects/galamsey-monotoring/assets/infn8vzn_mercury_training_340');
var roads = ee.FeatureCollection('projects/galamsey-monotoring/assets/ghana_osm_roads');

// Balanced training: match class 0 count to class 1 count
var class1 = trainingPoints.filter(ee.Filter.eq('class', 1));
var class1Sampled = imageStack.sampleRegions({
  collection: class1, properties: ['class'], scale: 10, tileScale: 8
});

var worldCover = ee.ImageCollection('ESA/WorldCover/v200').first().clip(AOI);

// Sample 500 class 0 points across all non-mining land cover types
var cleanPoints = imageStack
  .updateMask(worldCover.neq(50))   // exclude urban
  .updateMask(NDVI.gt(0.3))          // must have some vegetation
  .updateMask(BSI.lt(0.1))           // exclude bare soil
  .updateMask(MNDWI.lt(0.1))         // exclude water
  .sample({region: AOI, scale: 30, numPixels: 500, seed: 42, geometries: true})
  .map(function(f) { return f.set('class', 0); });

// Force equal class sizes — undersample class 0 to match class 1
var n = class1Sampled.size();
var cleanBalanced = cleanPoints.limit(n);

var allSamples = class1Sampled.merge(cleanBalanced).randomColumn('random', 42);
var trainSet = allSamples.filter(ee.Filter.lt('random', 0.8));
var valSet   = allSamples.filter(ee.Filter.gte('random', 0.8));

print('Class 1:', class1Sampled.size());
print('Class 0 balanced:', cleanBalanced.size());

var bands = ['B2','B3','B4','B5','B6','B7','B8','B8A','B11','B12',
             'NDVI','NDWI','MNDWI','NDTI','BSI','IronOxide'];

var classifier = ee.Classifier.smileRandomForest({numberOfTrees: 200, minLeafPopulation: 10, seed: 42})
  .train({features: trainSet, classProperty: 'class', inputProperties: bands})
  .setOutputMode('PROBABILITY');

var probability = imageStack.classify(classifier).rename('prob');

// Only flag pixels with >75% contamination probability
// NOTE: `finalMask` is referenced before its declaration below — works at GEE
// runtime due to function hoisting + late evaluation, but reads as a bug.
var contaminationOnly = probability.gt(0.75)
  .selfMask()
  .updateMask(finalMask)
  .rename('contamination');

Map.addLayer(probability, {min:0, max:1, palette:['white','orange','red']}, 'Contamination Probability', false);

// SHADOWING NOTE: the next two lines redeclare `classified` (= contaminationOnly above)
// and `worldCover` (already declared above). The second `contaminationOnly` reassignment
// also shadows the probability-thresholded version with the simpler classified.updateMask.
var classified = contaminationOnly;
var worldCover = ee.ImageCollection('ESA/WorldCover/v200').first().clip(AOI);
var urbanMask = worldCover.neq(50).and(worldCover.neq(40));
var roadBuffered = roads.map(function(f) { return f.buffer(50); });
var roadMask = ee.Image(1).byte().paint(roadBuffered, 0);
var finalMask = urbanMask.and(roadMask);
var contaminationOnly = classified.updateMask(finalMask);

var contaminationAreaDict = contaminationOnly
  .multiply(ee.Image.pixelArea()).divide(10000)
  .reduceRegion({reducer: ee.Reducer.sum(), geometry: AOI, scale: 100, maxPixels: 1e13, bestEffort: true});

// ── STEP 1: FLOOD SAFETY LAYER (JRC) ───────────────────────────────────
var jrc = ee.ImageCollection('JRC/GSW1_4/MonthlyHistory')
  .filter(ee.Filter.date('2000-01-01', '2021-12-31'));

var dry = jrc.filter(ee.Filter.calendarRange(11, 2, 'month'))
  .map(function(img) { return img.eq(2).selfMask(); })
  .max().clip(AOI);

var wet = jrc.filter(ee.Filter.calendarRange(6, 9, 'month'))
  .map(function(img) { return img.eq(2).selfMask(); })
  .max().clip(AOI);

var flood_zone = wet.unmask(0).subtract(dry.unmask(0))
  .gt(0).selfMask().clip(AOI);

var permanent_water = dry.unmask(0).and(wet.unmask(0)).selfMask().clip(AOI);

// ── BASIN GEOMETRIES ───────────────────────────────────────────────────
var basin_pra     = ee.Geometry.Polygon([[[-1.95,5.85],[-1.00,5.85],[-1.00,6.85],[-1.95,6.85],[-1.95,5.85]]], null, false);
var basin_ankobra = ee.Geometry.Polygon([[[-2.90,5.00],[-2.10,5.00],[-2.10,6.20],[-2.90,6.20],[-2.90,5.00]]], null, false);
var basin_birim   = ee.Geometry.Polygon([[[-1.45,5.65],[-0.45,5.65],[-0.45,6.45],[-1.45,6.45],[-1.45,5.65]]], null, false);
var basin_tano    = ee.Geometry.Polygon([[[-2.80,6.60],[-2.10,6.60],[-2.10,7.55],[-2.80,7.55],[-2.80,6.60]]], null, false);
var basin_offin   = ee.Geometry.Polygon([[[-2.20,5.95],[-1.55,5.95],[-1.55,6.70],[-2.20,6.70],[-2.20,5.95]]], null, false);
var basin_fc = ee.FeatureCollection([
  ee.Feature(basin_pra, {name:'Pra'}),
  ee.Feature(basin_ankobra, {name:'Ankobra'}),
  ee.Feature(basin_birim, {name:'Birim'}),
  ee.Feature(basin_tano, {name:'Tano'}),
  ee.Feature(basin_offin, {name:'Offin'})
]);

// ── FLOOD HA PER BASIN ─────────────────────────────────────────────────
var pixelArea = ee.Image.pixelArea();
var basins = [
  {name:'Pra', geom:basin_pra},
  {name:'Ankobra', geom:basin_ankobra},
  {name:'Birim', geom:basin_birim},
  {name:'Tano', geom:basin_tano},
  {name:'Offin', geom:basin_offin}
];

// ── VISUALIZE ALL LAYERS ───────────────────────────────────────────────
Map.addLayer(s2, {bands:['B4','B3','B2'], min:0, max:0.3}, 'True Colour', true);
Map.addLayer(permanent_water, {palette:['#003087']}, 'Permanent Water', true, 1.0);
Map.addLayer(dry, {palette:['#1a6eb5']}, 'Dry Season Water (Nov-Feb)', false);
Map.addLayer(wet, {palette:['#4da6ff']}, 'Wet Season Water (Jun-Sep)', false);
Map.addLayer(flood_zone, {palette:['#00ffff']}, 'Seasonal Flood Zone', true, 0.9);
Map.addLayer(turbidWater, {palette:['#FF6B00']}, 'Turbid Water (Mercury Indicator)', true, 0.85);
Map.addLayer(contaminationOnly, {palette:['#FF0000']}, 'Mercury Contamination (RF Classified)', true, 0.8);
Map.addLayer(trainingPoints, {color:'yellow'}, 'Training Points', false);
Map.addLayer(basin_fc.style({color:'FF6600', fillColor:'00000000', width:2}), {}, 'River Basins', true);
Map.addLayer(ee.Image().paint(ee.FeatureCollection([ee.Feature(AOI)]), 0, 2), {palette:['#FFFFFF']}, 'AOI', true);

// ── CONSOLE OUTPUT ────────────────────────────────────────────────────
print('=== GALASAT — ASHANTI REGION ===');
print('');
print('── TURBID WATER ──────────────────────');
print('Turbid water area (ha):', turbidAreaDict.get('area'));
print('');
print('── CONTAMINATION (RF) ────────────────');
print('Mercury contamination area (ha):', contaminationAreaDict.get('classification'));
print('');
print('── FLOOD SAFETY — BASIN EXTENT ───────');
basins.forEach(function(b) {
  var ha = ee.Number(
    pixelArea.updateMask(flood_zone.clip(b.geom))
      .reduceRegion({reducer:ee.Reducer.sum(), geometry:b.geom, scale:30, maxPixels:1e10, bestEffort:true})
      .get('area')
  ).divide(10000);
  print(ee.String(b.name).cat(': ').cat(ha.format('%.1f')).cat(' ha'));
});
print('SAFE: Nov-Feb | RISK: Jun-Sep');
print('');
print('── VALIDATION ────────────────────────');
var validated = valSet.classify(classifier);
var errorMatrix = validated.errorMatrix('class', 'classification');
print('Confusion Matrix:', errorMatrix);
print('Overall Accuracy:', errorMatrix.accuracy());
print('Kappa:', errorMatrix.kappa());
