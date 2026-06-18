// GalaSat v2.1 — INFN8 VZN + Community Search
// Mercury Contamination Intelligence
// Ashanti Region, Ghana | Remediation Sequencing + Biomonitoring Support
// Base: original working code — detector improved, SRTM added, collapsible legend
// Canonical v2.0 locked at: https://code.earthengine.google.com/6523f83ca8bfbf19a27431e16d50d11e
// DO NOT modify canonical — this is v2.1 working copy

// ============================================================
// STUDY AREA
// ============================================================
var trainingRef = ee.FeatureCollection(
  'projects/galamsey-monotoring/assets/infn8vzn_mercury_training_340'
);
var AOI = ee.Geometry.Rectangle([-2.60, 5.90, -1.00, 6.90]);
Map.setOptions('SATELLITE');
// Map.centerObject called after all layers added — faster initial tile load

// ============================================================
// POLYGON TRAINING ASSET
// ============================================================
var mappedPolygons = ee.FeatureCollection(
  'projects/galamsey-monotoring/assets/MINING_POLYGONS_800'
);

var labeledPolygons = mappedPolygons.map(function(f) {
  var name = ee.String(f.get('Name')).toUpperCase();
  var isMining = name.match('MERCURY').length().gt(0)
    .or(name.match('GALAMSEY').length().gt(0))
    .or(name.match('PIT').length().gt(0))
    .or(name.match('LAKE').length().gt(0));
  var isUrban = name.match('TOWN').length().gt(0)
    .or(name.match('HOUSE').length().gt(0))
    .or(name.match('VILLAGE').length().gt(0))
    .or(name.match('POPULATED').length().gt(0));
  var isScrub = name.match('SCRUB').length().gt(0)
    .or(name.match('DRY').length().gt(0))
    .or(name.match('BARE').length().gt(0))
    .or(name.match('DIRT').length().gt(0))
    .or(name.match('PATCH').length().gt(0))
    .or(name.match('FOREST').length().gt(0));
  var cls = ee.Algorithms.If(isMining, 1,
              ee.Algorithms.If(isUrban, 0,
                ee.Algorithms.If(isScrub, 2, -1)));
  return f.set('c', cls);
});

var miningPolys = labeledPolygons.filter(ee.Filter.eq('c', 1));
var urbanPolys  = labeledPolygons.filter(ee.Filter.eq('c', 0));
var scrubPolys  = labeledPolygons.filter(ee.Filter.eq('c', 2));

print('Mapped Polygon Classes');
print('Mining polygons:', miningPolys.size());
print('Urban polygons:', urbanPolys.size());
print('Scrubland polygons:', scrubPolys.size());

var miningPolyHa = ee.Number(miningPolys.map(function(f) {
  return f.set('area_ha', f.geometry().area(1).divide(10000));
}).aggregate_sum('area_ha'));
print('Manually mapped mining area (ha):', miningPolyHa.round());

// ============================================================
// OSM ASSETS
// ============================================================
var roadsAsset = ee.FeatureCollection(
  'projects/galamsey-monotoring/assets/ghana_osm_roads'
).filterBounds(AOI);
var waterways = ee.FeatureCollection(
  'projects/galamsey-monotoring/assets/ghana_osm_waterways'
).filterBounds(AOI);
var places = ee.FeatureCollection(
  'projects/galamsey-monotoring/assets/ghana_osm_places'
).filterBounds(AOI)
 .filter(ee.Filter.inList('fclass',[
   'village','hamlet','town','locality','suburb','isolated_dwelling','farm'
 ]));
var roads = roadsAsset.filter(ee.Filter.inList('fclass',
  ['motorway','trunk','primary','secondary'])).limit(500);

var tier1Rivers = waterways.filter(ee.Filter.inList('name',[
  'Pra','Pra River','Offin','Offin River','Ankobra','Ankobra River',
  'Birim','Birim River','Tano','Tano River','Densu','Densu River',
  'Oda','Oda River','Bonsa','Bonsa River','Subin','Subin River',
  'Afram','Afram River','Fena','Fena River','Anum','Anum River'
]));

// ============================================================
// GLOBAL MINING WATCH
// ============================================================
var globalMines = ee.FeatureCollection(
  'projects/sat-io/open-datasets/global-mining/global_mining_polygons'
).filterBounds(AOI);
var minesClipped = globalMines.map(function(f){
  return f.intersection(AOI, ee.ErrorMargin(1));
});
var totalMinesHa = ee.Number(minesClipped.map(function(f){
  return f.set('area_ha', f.geometry().area(1).divide(10000));
}).aggregate_sum('area_ha'));
var miningBufferImg = ee.Image(0).byte().paint(
  miningPolys.map(function(f){return f.buffer(500);}), 1
);

print('Mining Footprints');
print('GMW polygons in AOI:', globalMines.size());
print('GMW footprint area (ha):', totalMinesHa.round());
print('Confirmed river segments:', tier1Rivers.size());
print('Named places in AOI:', places.size());

// ============================================================
// CLOUD MASKING
// ============================================================
function maskS2(img) {
  var qa  = img.select('QA60');
  var qm  = qa.bitwiseAnd(1<<10).eq(0).and(qa.bitwiseAnd(1<<11).eq(0));
  var scl = img.select('SCL');
  var sm  = scl.neq(3).and(scl.neq(8)).and(scl.neq(9)).and(scl.neq(10)).and(scl.neq(11));
  return img.updateMask(qm).updateMask(sm).divide(10000);
}
function maskS2wet(img) {
  var scl = img.select('SCL');
  return img.updateMask(
    scl.neq(3).and(scl.neq(8)).and(scl.neq(9)).and(scl.neq(10))
  ).divide(10000);
}
function maskL8(img) {
  var qa = img.select('QA_PIXEL');
  return img.updateMask(
    qa.bitwiseAnd(1<<3).eq(0).and(qa.bitwiseAnd(1<<4).eq(0))
  ).multiply(0.0000275).add(-0.2);
}

// ============================================================
// COMPOSITES
// ============================================================
var s2dry = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
  .filterBounds(AOI).filterDate('2025-10-01','2026-02-28')
  .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE',30))
  .map(maskS2).median().clip(AOI);

var s2wet = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
  .filterBounds(AOI).filterDate('2023-05-01','2025-10-31')
  .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE',30))
  .map(maskS2wet).median().clip(AOI);

var now = ee.Date(Date.now());
var currentWeek = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
  .filterBounds(AOI)
  .filterDate(now.advance(-7,'day'), now)
  .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE',30))
  .map(maskS2).median().clip(AOI);

var previousWeek = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
  .filterBounds(AOI)
  .filterDate(now.advance(-14,'day'), now.advance(-7,'day'))
  .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE',30))
  .map(maskS2).median().clip(AOI);

var currentWeekSize = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
  .filterBounds(AOI).filterDate(now.advance(-7,'day'), now)
  .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE',30)).size();
var previousWeekSize = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
  .filterBounds(AOI).filterDate(now.advance(-14,'day'), now.advance(-7,'day'))
  .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE',30)).size();

print('Dry season images:',
  ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
    .filterBounds(AOI).filterDate('2025-10-01','2026-02-28')
    .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE',30)).size()
);
print('Current week images:', currentWeekSize);
print('Previous week images:', previousWeekSize);

// ============================================================
// SENTINEL-1 SAR
// ============================================================
function processSAR(image) {
  var vv = image.select('VV').focal_mean(3,'square','pixels').rename('VV_filtered');
  var vh = image.select('VH').focal_mean(3,'square','pixels').rename('VH_filtered');
  return vv.addBands(vh);
}

var s1dry = ee.ImageCollection('COPERNICUS/S1_GRD')
  .filterBounds(AOI).filterDate('2025-10-01','2026-03-31')
  .filter(ee.Filter.eq('instrumentMode','IW'))
  .filter(ee.Filter.listContains('transmitterReceiverPolarisation','VV'))
  .filter(ee.Filter.listContains('transmitterReceiverPolarisation','VH'))
  .filter(ee.Filter.eq('orbitProperties_pass','ASCENDING'))
  .select(['VV','VH']).map(processSAR).median().clip(AOI);

var s1wet = ee.ImageCollection('COPERNICUS/S1_GRD')
  .filterBounds(AOI).filterDate('2025-06-01','2025-09-30')
  .filter(ee.Filter.eq('instrumentMode','IW'))
  .filter(ee.Filter.listContains('transmitterReceiverPolarisation','VV'))
  .filter(ee.Filter.listContains('transmitterReceiverPolarisation','VH'))
  .filter(ee.Filter.eq('orbitProperties_pass','ASCENDING'))
  .select(['VV','VH']).map(processSAR).median().clip(AOI);

var sarRatio  = s1dry.select('VV_filtered').divide(s1dry.select('VH_filtered')).rename('SAR_ratio');
var sarChange = s1wet.select('VV_filtered').subtract(s1dry.select('VV_filtered')).rename('SAR_change');

print('SAR loaded - radar imagery active');

// ============================================================
// SPECTRAL INDICES
// ============================================================
function addIndicesS2(img) {
  var b2=img.select('B2'),b3=img.select('B3'),b4=img.select('B4');
  var b8=img.select('B8'),b8a=img.select('B8A'),b11=img.select('B11');
  return img
    .addBands(b11.add(b4).subtract(b8a.add(b2)).divide(b11.add(b4).add(b8a).add(b2)).rename('BSI'))
    .addBands(b8.subtract(b4).divide(b8.add(b4)).rename('NDVI'))
    .addBands(b3.subtract(b8).divide(b3.add(b8)).rename('NDWI'))
    .addBands(b11.subtract(b8).divide(b11.add(b8)).rename('NDBI'))
    .addBands(b3.subtract(b11).divide(b3.add(b11)).rename('MNDWI'))
    .addBands(b4.divide(b2).rename('IronOxide'));
}

// ============================================================
// COMPOSITE BUILDERS
// ============================================================
function compS2(d1,d2) {
  return ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
    .filterBounds(AOI).filterDate(d1,d2)
    .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE',30))
    .map(maskS2).median().clip(AOI);
}
function compL8(d1,d2) {
  return ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
    .filterBounds(AOI).filterDate(d1,d2).map(maskL8).median().clip(AOI);
}
function compL8rgb(d1,d2) {
  return compL8(d1,d2).select(['SR_B4','SR_B3','SR_B2'],['B4','B3','B2']);
}
function applyMinArea(mask,minPx) {
  return mask.updateMask(mask.connectedPixelCount(200,true).gte(minPx));
}

// ============================================================
// MASKS
// ============================================================
var worldCover = ee.ImageCollection('ESA/WorldCover/v200').mosaic().clip(AOI);
var urbanMaskBuffered = worldCover.neq(50).and(worldCover.neq(40))
  .focal_min(1,'square','pixels');
var ghsl = ee.ImageCollection('JRC/GHSL/P2016/SMOD_POP_GLOBE_V1')
  .filterBounds(AOI).mosaic().clip(AOI);
var ghslUrban = ghsl.select('smod_code').lt(3);
var roadMask = ee.Image(1).byte().paint(
  roads.map(function(f){return f.buffer(50);}), 0
);
var finalMask = urbanMaskBuffered.and(roadMask).and(ghslUrban);

// ============================================================
// RULE-BASED DETECTORS
// ============================================================
function mineS2(img) {
  var b2=img.select('B2'),b3=img.select('B3'),b4=img.select('B4');
  var b8=img.select('B8'),b8a=img.select('B8A'),b11=img.select('B11');
  var valid=b2.mask().and(b3.mask()).and(b4.mask()).and(b8.mask()).and(b8a.mask()).and(b11.mask());
  var bsi=b11.add(b4).subtract(b8a.add(b2)).divide(b11.add(b4).add(b8a).add(b2));
  var ndvi=b8.subtract(b4).divide(b8.add(b4));
  var ndwi=b3.subtract(b8).divide(b3.add(b8));
  var ndbi=b11.subtract(b8).divide(b11.add(b8));
  var mndwi=b3.subtract(b11).divide(b3.add(b11));
  var combined=
    bsi.gt(0.005).and(ndvi.lt(0.35)).and(ndbi.lt(0.05))
    .or(ndwi.gt(0.02).and(ndvi.lt(0.15)).and(bsi.gt(-0.08)))
    .or(b4.gt(0.07).and(b3.gt(0.05)).and(ndvi.lt(0.1)).and(ndbi.lt(0.0)))
    .or(mndwi.gt(0.0).and(ndvi.lt(0.1)).and(b4.gt(0.05)))
    .or(bsi.gt(0.05).and(ndvi.lt(0.45)))
    .or(b4.gt(0.12).and(ndvi.lt(0.2)).and(ndbi.lt(0.1)))
    .or(bsi.gt(0.02).and(ndbi.gt(0.05)).and(ndvi.lt(0.3)));
  return bsi.updateMask(applyMinArea(combined,6).and(valid).and(finalMask));
}

function mineL8(img) {
  var b2=img.select('SR_B2'),b3=img.select('SR_B3'),b4=img.select('SR_B4');
  var b5=img.select('SR_B5'),b6=img.select('SR_B6');
  var valid=b2.mask().and(b3.mask()).and(b4.mask()).and(b5.mask()).and(b6.mask());
  var bsi=b6.add(b4).subtract(b5.add(b2)).divide(b6.add(b4).add(b5).add(b2));
  var ndvi=b5.subtract(b4).divide(b5.add(b4));
  var ndwi=b3.subtract(b5).divide(b3.add(b5));
  var ndbi=b6.subtract(b5).divide(b6.add(b5));
  var mndwi=b3.subtract(b6).divide(b3.add(b6));
  var combined=bsi.gt(0.005).and(ndvi.lt(0.35)).and(ndbi.lt(0.05))
    .or(ndwi.gt(0.02).and(ndvi.lt(0.15)).and(bsi.gt(-0.08)))
    .or(b4.gt(0.07).and(b3.gt(0.05)).and(ndvi.lt(0.1)).and(ndbi.lt(0.0)))
    .or(mndwi.gt(0.0).and(ndvi.lt(0.1)).and(b4.gt(0.05)));
  return bsi.updateMask(applyMinArea(combined,6).and(valid).and(ghslUrban).and(roadMask));
}

function detectNewMining(current, previous) {
  var b2c=current.select('B2'),b3c=current.select('B3'),b4c=current.select('B4');
  var b8c=current.select('B8'),b8ac=current.select('B8A'),b11c=current.select('B11');
  var b2p=previous.select('B2'),b3p=previous.select('B3'),b4p=previous.select('B4');
  var b8p=previous.select('B8'),b8ap=previous.select('B8A'),b11p=previous.select('B11');
  var bsiC=b11c.add(b4c).subtract(b8ac.add(b2c)).divide(b11c.add(b4c).add(b8ac).add(b2c));
  var ndviC=b8c.subtract(b4c).divide(b8c.add(b4c));
  var bsiP=b11p.add(b4p).subtract(b8ap.add(b2p)).divide(b11p.add(b4p).add(b8ap).add(b2p));
  var ndviP=b8p.subtract(b4p).divide(b8p.add(b4p));
  var isMiningNow  = bsiC.gt(0.005).and(ndviC.lt(0.35)).and(finalMask);
  var wasNotBefore = bsiP.lt(0.005).or(ndviP.gt(0.35));
  return isMiningNow.and(wasNotBefore)
    .updateMask(isMiningNow.and(wasNotBefore).connectedPixelCount(200,true).gte(25))
    .selfMask().rename('NewMining');
}

// ============================================================
// RF FEATURE STACK
// ============================================================
var rfBands = ['B2','B3','B4','B5','B6','B7','B8','B8A','B11','B12',
               'BSI','NDVI','NDWI','NDBI','MNDWI','IronOxide',
               'VV_filtered','VH_filtered','SAR_ratio','SAR_change'];
var s2_rf      = compS2('2024-10-01','2025-04-30');
var s2_indexed = addIndicesS2(s2_rf);
var imageStack = s2_indexed.select(
  ['B2','B3','B4','B5','B6','B7','B8','B8A','B11','B12',
   'BSI','NDVI','NDWI','NDBI','MNDWI','IronOxide']
).addBands(s1dry.select('VV_filtered'))
 .addBands(s1dry.select('VH_filtered'))
 .addBands(sarRatio)
 .addBands(sarChange);

// IRON OXIDE ANOMALY
var ironOxideAnomalyThreshold = 2.5;
var ironOxideAnomaly = s2_indexed.select('IronOxide')
  .updateMask(s2_indexed.select('IronOxide').gt(ironOxideAnomalyThreshold))
  .updateMask(finalMask)
  .selfMask();

// ============================================================
// RF TRAINING
// ============================================================
var mSampled = imageStack.sampleRegions({
  collection:miningPolys, properties:['c'], scale:30, tileScale:4
}).randomColumn('r1',42).sort('r1').limit(1200);
var uSampled = imageStack.sampleRegions({
  collection:urbanPolys, properties:['c'], scale:30, tileScale:4
}).randomColumn('r2',42).sort('r2').limit(1200);
var sSampled = imageStack.sampleRegions({
  collection:scrubPolys, properties:['c'], scale:500, tileScale:4
});

var allSampled = mSampled.merge(uSampled).merge(sSampled);
var withRandom = allSampled.randomColumn('rand', 42);
var trainSet   = withRandom.filter(ee.Filter.lt('rand', 0.8));
var valSet     = withRandom.filter(ee.Filter.gte('rand', 0.8));

print('RF Training');
print('Mining samples:', mSampled.size());
print('Urban samples:', uSampled.size());
print('Scrubland samples:', sSampled.size());
print('Total samples:', allSampled.size());
print('Training set:', trainSet.size());
print('Validation set:', valSet.size());

var rfClassifier = ee.Classifier.smileRandomForest({numberOfTrees:200, seed:42})
  .train({features:trainSet, classProperty:'c', inputProperties:rfBands});
print('Classifier trained.');

var rfClass  = imageStack.classify(rfClassifier);
var ndviImg  = s2_indexed.select('NDVI');

var builtUp      = urbanMaskBuffered.not();
var denseUrban   = ghsl.select('smod_code').gt(2);
var urbanExclude = builtUp.or(denseUrban).unmask(0);
var roadsWithBuffer = roadsAsset.map(function(f) {
  var isMajor = ee.List(['motorway','trunk','primary']).contains(f.get('fclass'));
  var buf = ee.Number(ee.Algorithms.If(isMajor, 60, 25));
  return f.buffer(buf);
});
var roadFootprint = ee.Image(1).byte().paint(roadsWithBuffer, 0);

var rfRaw = rfClass.eq(1)
  .updateMask(rfClass.eq(1))
  .updateMask(urbanExclude.not())
  .updateMask(roadFootprint)
  .updateMask(ndviImg.lt(0.65));
var rfMining = rfRaw
  .updateMask(rfRaw.connectedPixelCount(200, true).gte(20))
  .selfMask();

// ============================================================
// VALIDATION
// ============================================================
var validated   = valSet.classify(rfClassifier);
var errorMatrix = validated.errorMatrix('c','classification');
print('RF Accuracy');
print('Confusion Matrix:', errorMatrix);
print('Overall Accuracy:', errorMatrix.accuracy());
print('Kappa:', errorMatrix.kappa());
print('Producers Accuracy:', errorMatrix.producersAccuracy());
print('Consumers Accuracy:', errorMatrix.consumersAccuracy());

// ============================================================
// TIME SERIES
// ============================================================
var l8_14=compL8('2013-10-01','2014-04-30');
var l8_15=compL8('2014-10-01','2015-04-30');
var s2_16=compS2('2015-06-01','2016-12-31');
var s2_18=compS2('2017-06-01','2018-06-30');
var s2_20=compS2('2019-10-01','2020-04-30');
var s2_22=compS2('2021-10-01','2022-04-30');
var s2_24=compS2('2023-10-01','2024-04-30');
var s2_25=compS2('2024-10-01','2025-04-30');
var s2_26=compS2('2025-10-01','2026-03-09');
var l8_16=compL8rgb('2015-06-01','2016-12-31');
var l8_18=compL8rgb('2017-06-01','2018-06-30');
var d16=s2_16.select(['B4','B3','B2']).unmask(l8_16).clip(AOI);
var d18=s2_18.select(['B4','B3','B2']).unmask(l8_18).clip(AOI);

var ml14=mineL8(l8_14),ml15=mineL8(l8_15);
var m16=mineS2(s2_16),m18=mineS2(s2_18),m20=mineS2(s2_20);
var m22=mineS2(s2_22),m24=mineS2(s2_24),m25=mineS2(s2_25),m26=mineS2(s2_26);

// ============================================================
// CUMULATIVE DISTURBANCE
// ============================================================
var cumulativeMining = ml14.unmask(0).gt(0)
  .or(ml15.unmask(0).gt(0)).or(m16.unmask(0).gt(0))
  .or(m18.unmask(0).gt(0)).or(m20.unmask(0).gt(0))
  .or(m22.unmask(0).gt(0)).or(m24.unmask(0).gt(0))
  .or(m25.unmask(0).gt(0)).or(m26.unmask(0).gt(0))
  .selfMask();

var newMiningActivity = ee.Algorithms.If(
  currentWeekSize.gt(0).and(previousWeekSize.gt(0)),
  detectNewMining(currentWeek, previousWeek),
  ee.Image(0).selfMask()
);
var newMiningImg = ee.Image(newMiningActivity).selfMask();
var newMiningHa = ee.Number(0);
print('New activity layer active — hectare calculation deferred to final build');

// ============================================================
// TURBID WATER
// ============================================================
var MNDWI_wet = s2wet.normalizedDifference(['B3','B11']).rename('MNDWI_wet');
var NDTI_wet  = s2wet.normalizedDifference(['B4','B3']).rename('NDTI_wet');
var IOR_wet   = s2wet.select('B4').divide(s2wet.select('B2')).rename('IOR_wet');
var NDVI_wet  = s2wet.normalizedDifference(['B8','B4']).rename('NDVI_wet');

var turbidWater = MNDWI_wet.gt(0.1).and(NDTI_wet.gt(0.08))
  .and(IOR_wet.gt(1.1)).and(NDVI_wet.lt(0.2))
  .and(s2wet.select('B3').gt(0.05))
  .selfMask().rename('TurbidWater');

// ============================================================
// FLOOD SAFETY JRC
// ============================================================
var pixelAreaImg  = ee.Image.pixelArea();
var jrcOccurrence = ee.Image('JRC/GSW1_4/GlobalSurfaceWater').select('occurrence').clip(AOI);
var jrcSeasonality= ee.Image('JRC/GSW1_4/GlobalSurfaceWater').select('seasonality').clip(AOI);
var permWater     = jrcOccurrence.gt(75).selfMask().clip(AOI);
var floodZone     = jrcSeasonality.gt(0).and(jrcOccurrence.lte(75)).selfMask().clip(AOI);
var forest        = worldCover.eq(10).selfMask();

// ============================================================
// BASIN GEOMETRIES
// ============================================================
var basin_pra     = ee.Geometry.Polygon([[[-1.95,5.85],[-1.00,5.85],[-1.00,6.85],[-1.95,6.85],[-1.95,5.85]]],null,false);
var basin_ankobra = ee.Geometry.Polygon([[[-2.90,5.00],[-2.10,5.00],[-2.10,6.20],[-2.90,6.20],[-2.90,5.00]]],null,false);
var basin_birim   = ee.Geometry.Polygon([[[-1.45,5.65],[-0.45,5.65],[-0.45,6.45],[-1.45,6.45],[-1.45,5.65]]],null,false);
var basin_tano    = ee.Geometry.Polygon([[[-3.20,6.50],[-2.00,6.50],[-2.00,7.80],[-3.20,7.80],[-3.20,6.50]]],null,false);
var basin_offin   = ee.Geometry.Polygon([[[-2.20,5.95],[-1.55,5.95],[-1.55,6.70],[-2.20,6.70],[-2.20,5.95]]],null,false);
var basin_fc = ee.FeatureCollection([
  ee.Feature(basin_pra,{name:'Pra'}),ee.Feature(basin_ankobra,{name:'Ankobra'}),
  ee.Feature(basin_birim,{name:'Birim'}),ee.Feature(basin_tano,{name:'Tano'}),
  ee.Feature(basin_offin,{name:'Offin'})
]);

// ============================================================
// HYDROSHEDS
// ============================================================
var hydroRivers       = ee.FeatureCollection('WWF/HydroSHEDS/v1/FreeFlowingRivers').filterBounds(AOI);
var hydroBasins       = ee.FeatureCollection('WWF/HydroSHEDS/v1/Basins/hybas_7').filterBounds(AOI);
var hydroRiversRaster = ee.Image().byte().paint({featureCollection:hydroRivers,color:1,width:2});
var hydroBasinsOutline= ee.Image().byte().paint({featureCollection:hydroBasins,color:1,width:2});

print('HydroSHEDS river segments:', hydroRivers.size());
print('HydroSHEDS watersheds:', hydroBasins.size());

// ============================================================
// SRTM ELEVATION + TERRAIN
// ============================================================
var srtm = ee.Image('USGS/SRTMGL1_003').clip(AOI);
var elevation = srtm.select('elevation');
var slope = ee.Terrain.slope(srtm);
var elevationSmooth = elevation.focal_mean(3,'square','pixels');
var downstreamZone = elevationSmooth.lt(150).and(slope.lt(5)).selfMask().rename('DownstreamZone');
var mercuryFlowRisk = downstreamZone
  .updateMask(miningBufferImg.unmask(0).focal_max(10,'square','pixels'))
  .selfMask().rename('MercuryFlowRisk');
var deploymentAccess = slope.lt(15).and(elevation.lt(300)).selfMask().rename('DeploymentAccess');

print('SRTM: elevation, slope, mercury flow risk, deployment access loaded');

// ============================================================
// LATERAL SPREAD MODEL
// ============================================================
var directSpreadZone = cumulativeMining.unmask(0).gt(0)
  .reproject({crs: 'EPSG:4326', scale: 100})
  .focal_max(5, 'circle', 'pixels')
  .reproject({crs: 'EPSG:4326', scale: 100})
  .selfMask()
  .rename('DirectSpread');

print('Lateral Spread Model — Layer A: Direct Spread — 500m buffer around confirmed mining');

var directSpreadHa = directSpreadZone.multiply(ee.Image.pixelArea()).divide(10000)
  .reduceRegion({
    reducer: ee.Reducer.sum(),
    geometry: AOI,
    scale: 100,
    maxPixels: 1e13,
    tileScale: 4,
    bestEffort: true
  });
print('Direct spread zone area (ha):', ee.Number(directSpreadHa.get('DirectSpread')).round());
print('Cumulative mining disturbance (ha): 170,960 — locked');

var riverProximityRisk = ee.Image(0).byte().paint(
  waterways.map(function(f){return f.buffer(100);}), 1
).selfMask();

var riverTransportRisk = riverProximityRisk.multiply(jrcSeasonality.divide(12).clamp(0,1))
  .add(riverProximityRisk.multiply(0.3))
  .clamp(0,1)
  .updateMask(miningBufferImg.unmask(0).focal_max(20,'square','pixels').gt(0))
  .selfMask()
  .rename('RiverTransportRisk');

var highRiverRisk = riverTransportRisk.gt(0.2).selfMask().rename('HighRiverRisk');

print('Lateral Spread Model — Layer B: River Transport Risk — downstream drainage network');

// ============================================================
// LANDSAT 5
// ============================================================
function maskL5(img) {
  var qa = img.select('QA_PIXEL');
  return img.updateMask(
    qa.bitwiseAnd(1<<3).eq(0).and(qa.bitwiseAnd(1<<4).eq(0))
  ).multiply(0.0000275).add(-0.2);
}
function compL5(d1,d2) {
  return ee.ImageCollection('LANDSAT/LT05/C02/T1_L2')
    .filterBounds(AOI).filterDate(d1,d2)
    .filter(ee.Filter.lt('CLOUD_COVER',30))
    .map(maskL5).median().clip(AOI);
}
function mineL5(img) {
  var bandCount = img.bandNames().size();
  var safeResult = ee.Algorithms.If(
    bandCount.gt(0),
    ee.Image(function() {
      var b1=img.select('SR_B1'),b2=img.select('SR_B2'),b3=img.select('SR_B3');
      var b4=img.select('SR_B4'),b5=img.select('SR_B5'),b7=img.select('SR_B7');
      var valid=b1.mask().and(b2.mask()).and(b3.mask()).and(b4.mask()).and(b5.mask()).and(b7.mask());
      var bsi=b7.add(b3).subtract(b4.add(b1)).divide(b7.add(b3).add(b4).add(b1));
      var ndvi=b4.subtract(b3).divide(b4.add(b3));
      var ndwi=b2.subtract(b4).divide(b2.add(b4));
      var ndbi=b5.subtract(b4).divide(b5.add(b4));
      var mndwi=b2.subtract(b5).divide(b2.add(b5));
      var combined=bsi.gt(0.005).and(ndvi.lt(0.35)).and(ndbi.lt(0.05))
        .or(ndwi.gt(0.02).and(ndvi.lt(0.15)).and(bsi.gt(-0.08)))
        .or(b3.gt(0.07).and(b2.gt(0.05)).and(ndvi.lt(0.1)).and(ndbi.lt(0.0)))
        .or(mndwi.gt(0.0).and(ndvi.lt(0.1)).and(b3.gt(0.05)));
      return bsi.updateMask(applyMinArea(combined,6).and(valid).and(ghslUrban).and(roadMask));
    }()),
    ee.Image(0).selfMask()
  );
  return ee.Image(safeResult);
}

print('L5 1984 images:', ee.ImageCollection('LANDSAT/LT05/C02/T1_L2').filterBounds(AOI).filterDate('1984-10-01','1985-04-30').size());
print('L5 1988 images:', ee.ImageCollection('LANDSAT/LT05/C02/T1_L2').filterBounds(AOI).filterDate('1987-10-01','1988-04-30').size());

var l5_84=compL5('1984-10-01','1985-04-30');
var l5_88=compL5('1987-10-01','1988-04-30');
var m84=mineL5(l5_84);
var m88=mineL5(l5_88);

var cumulativeMiningExtended = m84.unmask(0).gt(0)
  .or(m88.unmask(0).gt(0))
  .or(cumulativeMining.unmask(0).gt(0))
  .selfMask();

var extCumAreaDict = cumulativeMiningExtended.multiply(ee.Image.pixelArea()).divide(10000)
  .reduceRegion({reducer:ee.Reducer.sum(),geometry:AOI,scale:200,maxPixels:1e13,bestEffort:true});

print('Landsat 5 Extension 1984-2012');
print('Extended cumulative disturbance 1984-2026:',
  ee.Number(extCumAreaDict.values().get(0)).round());

// ============================================================
// WORLDPOP
// ============================================================
var worldPop = ee.ImageCollection('WorldPop/GP/100m/pop')
  .filterBounds(AOI)
  .filter(ee.Filter.eq('year', 2020))
  .mosaic()
  .clip(AOI);

var miningBuffer2km = cumulativeMining.unmask(0)
  .focal_max(20, 'square', 'pixels')
  .gt(0);

var exposedPop = worldPop.updateMask(miningBuffer2km);

var exposedPopTotal = exposedPop
  .reduceRegion({
    reducer: ee.Reducer.sum(),
    geometry: AOI,
    scale: 100,
    maxPixels: 1e13,
    tileScale: 4,
    bestEffort: true
  });

print('WorldPop Population Analysis');

// ============================================================
// GPM PRECIPITATION
// ============================================================
var gpmWetSeason = ee.ImageCollection('NASA/GPM_L3/IMERG_V06')
  .filterBounds(AOI)
  .filter(ee.Filter.calendarRange(6, 9, 'month'))
  .filterDate('2023-01-01','2025-12-31')
  .select('precipitationCal')
  .mean()
  .clip(AOI);

var gpmDrySeason = ee.ImageCollection('NASA/GPM_L3/IMERG_V06')
  .filterBounds(AOI)
  .filter(ee.Filter.calendarRange(11, 2, 'month'))
  .filterDate('2023-01-01','2025-12-31')
  .select('precipitationCal')
  .mean()
  .clip(AOI);

print('GPM Precipitation loaded — wet and dry season composites ready');

// ============================================================
// VIIRS NIGHT LIGHTS
// ============================================================
var viirs = ee.ImageCollection('NOAA/VIIRS/DNB/MONTHLY_V1/VCMSLCFG')
  .filterBounds(AOI)
  .filterDate('2024-01-01','2024-12-31')
  .select('avg_rad')
  .mean()
  .clip(AOI);

print('VIIRS Night Lights loaded — 2024 annual composite');
print('Wet season (Jun-Sep): mean precipitation mm/hr over AOI');
print('Dry season (Nov-Feb): mean precipitation mm/hr over AOI');
print('Total population within 2km of mining disturbance:',
  ee.Number(exposedPopTotal.values().get(0)).round());

var highDensityExposed = exposedPop.gt(50).selfMask();

// ============================================================
// BOREHOLE EXCLUSION ZONE
// ============================================================
var globgm = ee.Image('projects/sat-io/open-datasets/GLOBGM/STEADY-STATE/globgm-wtd-ss').clip(AOI);
var shallowWaterTable = globgm.lt(5).selfMask();
var boreholeExclusionZone = shallowWaterTable
  .updateMask(
    miningBufferImg.unmask(0).focal_max(30,'square','pixels')
    .or(turbidWater.unmask(0))
    .or(ee.Image(0).byte().paint(tier1Rivers.map(function(f){return f.buffer(1000);}),1).unmask(0))
  )
  .selfMask()
  .rename('BoreholeExclusion');

print('Borehole Exclusion Zone — shallow water table + contamination proximity');

// ============================================================
// FARMING SAFETY LAYER
// ============================================================
var croplandMask = worldCover.eq(40);
var contaminationEnvelope = miningBufferImg.unmask(0)
  .or(turbidWater.unmask(0))
  .or(ee.Image(0).byte().paint(
    tier1Rivers.map(function(f){return f.buffer(500);}), 1
  ).unmask(0));

var cropAtRisk = croplandMask.updateMask(contaminationEnvelope).selfMask().rename('CropAtRisk');
var cropSafe   = croplandMask.updateMask(contaminationEnvelope.not()).selfMask().rename('CropSafe');

var cropAtRiskHa = cropAtRisk.multiply(ee.Image.pixelArea()).divide(10000)
  .reduceRegion({
    reducer: ee.Reducer.sum(),
    geometry: AOI,
    scale: 500,
    maxPixels: 1e13,
    tileScale: 4,
    bestEffort: true
  });
print('Farming Safety — Cropland at mercury contamination risk (ha):',
  ee.Number(cropAtRiskHa.get('CropAtRisk')).round());
print('Farming Safety Layer: at-risk cropland (pink) | safe cropland (white)');

// ============================================================
// WATER RASTER LAYERS
// ============================================================
var allWaterwaysRaster      = ee.Image().byte().paint({featureCollection:waterways,color:1,width:1});
var tier1Raster             = ee.Image().byte().paint({featureCollection:tier1Rivers,color:1,width:3});
var adjacentWaterwaysRaster = ee.Image().byte().paint(
  waterways.map(function(f){return f.buffer(30);}),1
).updateMask(miningBufferImg);

// ============================================================
// COMMUNITY EXPOSURE
// ============================================================
var gmwMiningUnion = globalMines.map(function(f){return f.buffer(2000);})
  .union(ee.ErrorMargin(100));
var exposedCommunitiesGMW = places.filterBounds(gmwMiningUnion.geometry());

var miningFootprintUnion = globalMines.merge(miningPolys)
  .map(function(f){return f.buffer(2000, ee.ErrorMargin(10));})
  .union(ee.ErrorMargin(10));

var placesWithBuffer = places.filterBounds(miningFootprintUnion.geometry());

var exposedCommunities = ee.FeatureCollection(placesWithBuffer.map(function(f) {
  var pt = f.geometry();
  var basinName = ee.String(ee.Algorithms.If(
    basin_fc.filterBounds(pt).size().gt(0),
    basin_fc.filterBounds(pt).first().get('name'),
    'Ashanti'
  ));
  return f.set('basin_name', basinName);
}));

print('Community Exposure — GalaSat vs GMW buffer comparison');
print('Communities within 2km of GMW polygons:', exposedCommunitiesGMW.size());
print('Communities within 2km of GalaSat detection:', exposedCommunities.size());

// ============================================================
// COMMUNITY RISK SCORE
// ============================================================
var w_turbid  = 0.30;
var w_river   = 0.25;
var w_mining  = 0.25;
var w_flood   = 0.15;
var w_pop     = 0.05;

var turbidRiskImg = turbidWater.unmask(0).clamp(0,1).rename('turbid_risk');
var riverRiskImg = ee.Image(0).byte().paint(
  tier1Rivers.map(function(f){return f.buffer(500);}), 1
).clamp(0,1).rename('river_risk');
var miningRiskImg = miningBufferImg.unmask(0)
  .focal_max(10,'square','pixels').clamp(0,1).rename('mining_risk');
var floodRiskImg = jrcSeasonality.divide(12).clamp(0,1).rename('flood_risk');
var popRiskImg = worldPop.divide(200).clamp(0,1).rename('pop_risk');

var compositeRiskImg = turbidRiskImg.multiply(w_turbid)
  .add(riverRiskImg.multiply(w_river))
  .add(miningRiskImg.multiply(w_mining))
  .add(floodRiskImg.multiply(w_flood))
  .add(popRiskImg.multiply(w_pop))
  .clamp(0,1)
  .rename('risk_score');

print('Community Risk Model — Placeholder Weights (pending Pure Earth calibration)');
print('Turbid water: ' + (w_turbid*100) + '% | River proximity: ' + (w_river*100) + '% | Mining: ' + (w_mining*100) + '% | Flood: ' + (w_flood*100) + '% | Population: ' + (w_pop*100) + '%');

var bufferedCommunities = exposedCommunities.map(function(f) {
  return f.buffer(500);
});

var scoredCommunities = compositeRiskImg.addBands(miningRiskImg)
  .addBands(turbidRiskImg)
  .addBands(riverRiskImg)
  .addBands(popRiskImg)
  .addBands(floodRiskImg)
  .addBands(deploymentAccess.unmask(0).rename('deploy_access'))
  .reduceRegions({
    collection: bufferedCommunities,
    reducer: ee.Reducer.mean(),
    scale: 100,
    tileScale: 8,
    crs: 'EPSG:4326'
  });

var rankedCommunities = scoredCommunities.sort('risk_score', false);
var top100 = rankedCommunities.limit(100);
var top10  = rankedCommunities.limit(10);

print('Community Risk Score — Top 10 Highest Risk:');
print(top10);

var riskViz = {min:0,max:1,palette:['ffffff','b3e5fc','4fc3f7','0288d1','01579b']};
var communityRiskImg = rankedCommunities.reduceToImage({
  properties:['risk_score'],
  reducer:ee.Reducer.first()
});

// ============================================================
// AREA CALCULATIONS
// ============================================================
var turbidAreaDict = turbidWater.multiply(pixelAreaImg).divide(10000)
  .reduceRegion({reducer:ee.Reducer.sum(),geometry:AOI,scale:20,maxPixels:1e13,bestEffort:true});
print('Turbid Water Area (ha):', ee.Number(turbidAreaDict.get('TurbidWater')).round());

var rfAreaDict = rfMining.multiply(pixelAreaImg).divide(10000)
  .reduceRegion({reducer:ee.Reducer.sum(),geometry:AOI,scale:200,maxPixels:1e13,bestEffort:true});
print('RF Mining Ha:', ee.Number(rfAreaDict.get(rfAreaDict.keys().get(0))).round());

var cumAreaDict = cumulativeMining.multiply(pixelAreaImg).divide(10000)
  .reduceRegion({reducer:ee.Reducer.sum(),geometry:AOI,scale:200,maxPixels:1e13,bestEffort:true});
print('Cumulative Mining Ha (2014-2026):', ee.Number(cumAreaDict.get(cumAreaDict.keys().get(0))).round());

var mercuryFlowRiskHa = mercuryFlowRisk.multiply(pixelAreaImg).divide(10000)
  .reduceRegion({reducer:ee.Reducer.sum(),geometry:AOI,scale:30,maxPixels:1e13,bestEffort:true});
print('Mercury Flow Risk (ha):', ee.Number(mercuryFlowRiskHa.get('MercuryFlowRisk')).round());

print('Communities in AOI:', places.size());
print('Communities within 2km of mining:', exposedCommunities.size());

print('Flood Safety per basin:');
var basins=[
  {name:'Pra',geom:basin_pra},{name:'Ankobra',geom:basin_ankobra},
  {name:'Birim',geom:basin_birim},{name:'Tano',geom:basin_tano},
  {name:'Offin',geom:basin_offin}
];
basins.forEach(function(b){
  var ha=ee.Number(pixelAreaImg.updateMask(floodZone.clip(b.geom))
    .reduceRegion({reducer:ee.Reducer.sum(),geometry:b.geom,scale:30,maxPixels:1e10,bestEffort:true})
    .get('area')).divide(10000);
  print(ee.String(b.name).cat(': ').cat(ha.format('%.1f')).cat(' ha'));
});
print('SAFE deployment: Nov-Feb | RISK: Jun-Sep');

// ============================================================
// TIME SERIES CHART
// ============================================================
function getYearTurbidHa(year) {
  var s=year+'-06-01', e=year+'-09-30';
  function mL8(img){
    var qa=img.select('QA_PIXEL');
    return img.updateMask(qa.bitwiseAnd(1<<3).eq(0).and(qa.bitwiseAnd(1<<4).eq(0)))
      .select(['SR_B3','SR_B4','SR_B6']).multiply(0.0000275).add(-0.2);
  }
  var merged=ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
    .filterBounds(AOI).filterDate(s,e).filter(ee.Filter.lt('CLOUD_COVER',75)).map(mL8)
    .merge(ee.ImageCollection('LANDSAT/LC09/C02/T1_L2')
      .filterBounds(AOI).filterDate(s,e).filter(ee.Filter.lt('CLOUD_COVER',75)).map(mL8));
  var s2y=ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
    .filterBounds(AOI).filterDate(s,e).filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE',75)).map(maskS2wet);
  var s2Ha=ee.Algorithms.If(s2y.size().gte(3),
    ee.Number(s2y.median().clip(AOI).normalizedDifference(['B3','B11']).gt(0)
      .and(s2y.median().clip(AOI).normalizedDifference(['B4','B3']).gt(0.05))
      .multiply(pixelAreaImg).divide(10000)
      .reduceRegion({reducer:ee.Reducer.sum(),geometry:AOI,scale:30,maxPixels:1e13,bestEffort:true})
      .values().get(0)),0);
  var lHa=ee.Algorithms.If(merged.size().gte(3),
    ee.Number(merged.median().clip(AOI).normalizedDifference(['SR_B3','SR_B6']).gt(0)
      .and(merged.median().clip(AOI).normalizedDifference(['SR_B4','SR_B3']).gt(0.05))
      .multiply(pixelAreaImg).divide(10000)
      .reduceRegion({reducer:ee.Reducer.sum(),geometry:AOI,scale:30,maxPixels:1e13,bestEffort:true})
      .values().get(0)),null);
  return ee.Algorithms.If(merged.size().gte(3),lHa,s2Ha);
}

var tsFC=ee.FeatureCollection(
  [2014,2015,2016,2017,2018,2019,2020,2021,2022,2023,2024,2025].map(function(y){
    return ee.Feature(null,{year:y,turbid_ha:getYearTurbidHa(y)});
  })
);

print(ui.Chart.feature.byFeature({features:tsFC,xProperty:'year',yProperties:['turbid_ha']})
  .setChartType('LineChart').setOptions({
    title:'GalaSat v2.1 — INFN8 VZN | Mercury-Affected Turbid Water 2014–2025 | Ashanti Region, Ghana',
    hAxis:{title:'Year',format:'####'},
    vAxis:{title:'Turbid Water Area (hectares)'},
    colors:['#FF6D00'],lineWidth:3,pointSize:6,interpolateNulls:false
  }));

// ============================================================
// FOUR-TIER COMMUNITY DOTS
// ============================================================
var tier1Dots = ee.Image().byte().paint(
  rankedCommunities.limit(10).map(function(f){return f.buffer(800);}), 1
).selfMask();
var tier2Dots = ee.Image().byte().paint(
  rankedCommunities.limit(50).map(function(f){return f.buffer(600);}), 1
).selfMask();
var tier3Dots = ee.Image().byte().paint(
  rankedCommunities.limit(100).map(function(f){return f.buffer(400);}), 1
).selfMask();
var tier4Dots = ee.Image().byte().paint(
  exposedCommunities.map(function(f){return f.buffer(500);}), 1
).selfMask();

// ============================================================
// MAP LAYERS — store references for toggle tracking
// ============================================================
var layerRefs = [];

// AOI BOUNDARY — added first so it renders immediately on load
var aoiOutline = ee.Image().byte().paint({
  featureCollection: ee.FeatureCollection([ee.Feature(AOI)]), color:1, width:2
});
layerRefs.push(Map.addLayer(aoiOutline,{palette:['FFFFFF'],opacity:1.0},'Pilot AOI Boundary'));

// GROUP 0 — HISTORICAL RECORD (bottom of panel)
layerRefs.push(Map.addLayer(l5_84,{bands:['SR_B3','SR_B2','SR_B1'],min:0.0,max:0.3,gamma:1.4},'Historical — 1984 True Color (L5)',false));
layerRefs.push(Map.addLayer(m84,{min:0.0,max:0.3,palette:['FF6D00','FF8C00','FFA726','FFB74D']},'Historical — 1984 Mining (L5)',false));
layerRefs.push(Map.addLayer(l5_88,{bands:['SR_B3','SR_B2','SR_B1'],min:0.0,max:0.3,gamma:1.4},'Historical — 1988 True Color (L5)',false));
layerRefs.push(Map.addLayer(m88,{min:0.0,max:0.3,palette:['FF6D00','FF8C00','FFA726','FFB74D']},'Historical — 1988 Mining (L5)',false));
layerRefs.push(Map.addLayer(l8_14,{bands:['SR_B4','SR_B3','SR_B2'],min:0.0,max:0.3,gamma:1.4},'Historical — 2014 True Color',false));
layerRefs.push(Map.addLayer(ml14,{min:0.0,max:0.3,palette:['FF6D00','FF8C00','FFA726','FFB74D']},'Historical — 2014 Mining (L8)',false));
layerRefs.push(Map.addLayer(l8_15,{bands:['SR_B4','SR_B3','SR_B2'],min:0.0,max:0.3,gamma:1.4},'Historical — 2015 True Color',false));
layerRefs.push(Map.addLayer(ml15,{min:0.0,max:0.3,palette:['FF6D00','FF8C00','FFA726','FFB74D']},'Historical — 2015 Mining (L8)',false));
layerRefs.push(Map.addLayer(d16,{bands:['B4','B3','B2'],min:0.0,max:0.3,gamma:1.4},'Historical — 2016 True Color',false));
layerRefs.push(Map.addLayer(m16,{min:0.0,max:0.3,palette:['FF6D00','FF8C00','FFA726','FFB74D']},'Historical — 2016 Mining',false));
layerRefs.push(Map.addLayer(d18,{bands:['B4','B3','B2'],min:0.0,max:0.3,gamma:1.4},'Historical — 2018 True Color',false));
layerRefs.push(Map.addLayer(m18,{min:0.0,max:0.3,palette:['FF6D00','FF8C00','FFA726','FFB74D']},'Historical — 2018 Mining',false));
layerRefs.push(Map.addLayer(s2_20,{bands:['B4','B3','B2'],min:0.0,max:0.3,gamma:1.4},'Historical — 2020 True Color',false));
layerRefs.push(Map.addLayer(m20,{min:0.0,max:0.3,palette:['FF6D00','FF8C00','FFA726','FFB74D']},'Historical — 2020 Mining',false));
layerRefs.push(Map.addLayer(s2_22,{bands:['B4','B3','B2'],min:0.0,max:0.3,gamma:1.4},'Historical — 2022 True Color',false));
layerRefs.push(Map.addLayer(m22,{min:0.0,max:0.3,palette:['FF6D00','FF8C00','FFA726','FFB74D']},'Historical — 2022 Mining',false));
layerRefs.push(Map.addLayer(s2_24,{bands:['B4','B3','B2'],min:0.0,max:0.3,gamma:1.4},'Historical — 2024 True Color',false));
layerRefs.push(Map.addLayer(m24,{min:0.0,max:0.3,palette:['FF6D00','FF8C00','FFA726','FFB74D']},'Historical — 2024 Mining',false));
layerRefs.push(Map.addLayer(s2_25,{bands:['B4','B3','B2'],min:0.0,max:0.3,gamma:1.4},'Historical — 2025 True Color',false));
layerRefs.push(Map.addLayer(m25,{min:0.0,max:0.3,palette:['FF6D00','FF8C00','FFA726','FFB74D']},'Historical — 2025 Mining',false));
layerRefs.push(Map.addLayer(s2_26,{bands:['B4','B3','B2'],min:0.0,max:0.3,gamma:1.4},'Historical — 2026 True Color',false));
layerRefs.push(Map.addLayer(cumulativeMiningExtended,{min:0,max:1,palette:['90A4AE'],opacity:0.7},'Historical — Cumulative Mining Extended 1984-2026',false));

// GROUP 7 — REFERENCE
layerRefs.push(Map.addLayer(basin_fc.style({color:'B0BEC5',fillColor:'00000000',width:2}),{},'Reference — Basin Boundaries',false));
layerRefs.push(Map.addLayer(hydroBasinsOutline,{palette:['FFD600']},'Reference — HydroSHEDS Watersheds',false));
layerRefs.push(Map.addLayer(hydroRiversRaster,{palette:['0091EA']},'Reference — HydroSHEDS Rivers',false));
layerRefs.push(Map.addLayer(roads,{color:'616161'},'Reference — Major Roads',false));
layerRefs.push(Map.addLayer(sarChange,{min:-5,max:5,palette:['0000ff','ffffff','ff0000']},'Reference — SAR Change Wet-Dry',false));
layerRefs.push(Map.addLayer(s1dry.select('VV_filtered'),{min:-20,max:0,palette:['000000','ffffff']},'Reference — SAR VV Dry Season',false));
layerRefs.push(Map.addLayer(s2dry,{bands:['B4','B3','B2'],min:0.0,max:0.3,gamma:1.4},'Reference — Base Image Nov2025-Feb2026',false));

// GROUP 6 — ENVIRONMENT
layerRefs.push(Map.addLayer(cropAtRisk,{palette:['FF1493']},'Environment — Cropland at Mercury Risk (food chain)',false));
layerRefs.push(Map.addLayer(cropSafe,{palette:['FFFFFF']},'Environment — Cropland Safe from Contamination',false));
layerRefs.push(Map.addLayer(downstreamZone,{palette:['00897B'],opacity:0.6},'Environment — Downstream Low-Elevation Zones',false));
layerRefs.push(Map.addLayer(elevation,{min:0,max:400,palette:['1B5E20','4CAF50','FFEB3B','FF5722','B71C1C']},'Environment — Elevation SRTM 30m',false));
layerRefs.push(Map.addLayer(forest,{palette:['1B5E20']},'Environment — Forest Cover',false));
layerRefs.push(Map.addLayer(gpmDrySeason,{min:0,max:0.1,palette:['ffffff','ffffd4','fed98e','fe9929','d95f0e','993404']},'Environment — Precipitation Dry Season Nov-Feb (GPM)',false));
layerRefs.push(Map.addLayer(gpmWetSeason,{min:0,max:0.5,palette:['ffffff','e0f3f8','abd9e9','74add1','4575b4','313695']},'Environment — Precipitation Wet Season Jun-Sep (GPM)',false));
layerRefs.push(Map.addLayer(slope,{min:0,max:30,palette:['ffffff','ffcc00','ff6600','cc0000']},'Environment — Slope Terrain Grade',false));

// GROUP 5 — DEPLOYMENT & ACCESS
layerRefs.push(Map.addLayer(boreholeExclusionZone,{palette:['FF3D00'],opacity:0.75},'Deployment — Borehole Exclusion Zone (contamination risk)',false));
layerRefs.push(Map.addLayer(floodZone,{palette:['00BCD4']},'Deployment — Flood Zone Seasonal Safety',false,0.85));
layerRefs.push(Map.addLayer(deploymentAccess,{palette:['76FF03'],opacity:0.5},'Deployment — Remediation-Accessible Terrain',false));

// GROUP 4 — COMMUNITY EXPOSURE
layerRefs.push(Map.addLayer(viirs,{min:0,max:10,palette:['000000','1a1a2e','FFF9C4','FFEE58','FFD600']},'Community — Economic Activity Night Lights 2024 (VIIRS)',false));
layerRefs.push(Map.addLayer(highDensityExposed,{palette:['FF6F00'],opacity:0.8},'Community — High-Density Exposure Zones (>50 per pixel)',false));
layerRefs.push(Map.addLayer(worldPop,{min:0,max:200,palette:['ffffff','ffeda0','feb24c','f03b20','bd0026']},'Community — Population Density WorldPop 2020',false));
layerRefs.push(Map.addLayer(exposedPop,{min:0,max:200,palette:['fff3e0','FF6F00','e65100']},'Community — Population within 2km of Mining',false));
layerRefs.push(Map.addLayer(communityRiskImg,riskViz,'Community — Mercury Exposure Risk Score',false));
layerRefs.push(Map.addLayer(tier4Dots,{palette:['76FF03'],opacity:0.8},'Community — All 2km Communities Unranked',true));
layerRefs.push(Map.addLayer(tier3Dots,{palette:['0D47A1'],opacity:0.9},'Community — Biomonitoring Elevated Risk Ranks 51-100',true));
layerRefs.push(Map.addLayer(tier2Dots,{palette:['00E5FF'],opacity:0.9},'Community — Biomonitoring High Risk Ranks 11-50',true));
layerRefs.push(Map.addLayer(tier1Dots,{palette:['FF00FF'],opacity:1.0},'Community — Biomonitoring Critical Risk Ranks 1-10',true));

// GROUP 3 — MERCURY TRANSPORT
layerRefs.push(Map.addLayer(directSpreadZone,{palette:['D500F9'],opacity:0.7},'Mercury Transport — Direct Spread Zone 500m from Mining',false));
layerRefs.push(Map.addLayer(mercuryFlowRisk,{palette:['F9A825'],opacity:0.75},'Mercury Transport — Pathways SRTM-derived',false));
layerRefs.push(Map.addLayer(highRiverRisk,{palette:['7B1FA2'],opacity:0.75},'Mercury Transport — River Transport Risk Downstream Network',false));

// GROUP 2 — WATER CONTAMINATION
layerRefs.push(Map.addLayer(allWaterwaysRaster,{palette:['651FFF'],opacity:0.4},'Water — All Waterways',false));
layerRefs.push(Map.addLayer(tier1Raster,{palette:['8B0000']},'Water — Contaminated Rivers Field-Verified',false));
layerRefs.push(Map.addLayer(permWater,{palette:['FFFFFF']},'Water — Uncontaminated Water Bodies',false));
layerRefs.push(Map.addLayer(adjacentWaterwaysRaster,{palette:['651FFF']},'Water — Waterways Adjacent to Mining',false));
layerRefs.push(Map.addLayer(turbidWater,{palette:['8B0000'],opacity:0.8},'Water — Mercury-Contaminated Turbid Wet Season',false));

// GROUP 1 — MINING DISTURBANCE (top of panel)
layerRefs.push(Map.addLayer(urbanPolys.style({color:'FFC107',fillColor:'FFC10740',width:1}),{},'Mining — Mapped Polygons Urban',false));
layerRefs.push(Map.addLayer(scrubPolys.style({color:'00ff00',fillColor:'00ff0040',width:1}),{},'Mining — Mapped Polygons Scrubland',false));
layerRefs.push(Map.addLayer(miningPolys.style({color:'FF0000',fillColor:'FF000040',width:1}),{},'Mining — Mapped Polygons Active Sites',false));
layerRefs.push(Map.addLayer(globalMines,{color:'FFFF00'},'Mining — Footprints Global Mining Watch',false));
layerRefs.push(Map.addLayer(ironOxideAnomaly,{min:2.5,max:5.0,palette:['ff6f00','e53935','b71c1c','4a0000']},'Mining — Iron Oxide Anomaly ASGM Soil Signal',false));
layerRefs.push(Map.addLayer(rfMining,{min:0,max:1,palette:['7B1FA2'],opacity:0.7},'Mining — AI-Detected Disturbance 2025',false));
layerRefs.push(Map.addLayer(newMiningImg,{min:0,max:1,palette:['FF0000'],opacity:0.9},'Mining — New Activity This Week',true));
layerRefs.push(Map.addLayer(cumulativeMining,{min:0,max:1,palette:['FF6D00'],opacity:0.7},'Mining — Cumulative Disturbance 2014-2026',false));
layerRefs.push(Map.addLayer(m26,{min:0.0,max:0.3,palette:['FF6D00','FF8C00','FFA726','FFB74D']},'Mining — Current Activity Sentinel-2 Weekly Composite',true));

// Center map after all layers registered — faster tile load
Map.centerObject(AOI, 10);

// AOI boundary added first — see top of layer section

// ============================================================
// BIOMONITORING PRIORITY PANEL
// ============================================================
var riskPanelVisible = false;
var riskContent = ui.Panel({style:{margin:'0',padding:'0',shown:false}});
riskContent.add(ui.Label({value:'GalaSat biomonitoring prioritization model',style:{fontSize:'11px',color:'#333333',fontWeight:'bold',margin:'0 0 2px 0'}}));
riskContent.add(ui.Label({value:'Weights pending calibration against Pure Earth field mercury data',style:{fontSize:'9px',color:'#8B0000',margin:'0 0 2px 0',fontStyle:'italic'}}));
riskContent.add(ui.Label({value:'Turbid water 30% | River proximity 25% | Mining 25% | Flood 15% | Population 5%',style:{fontSize:'9px',color:'#888888',margin:'0 0 6px 0'}}));
var loadingLabel = ui.Label({value:'Loading community risk scores...',style:{fontSize:'11px',color:'#555555',margin:'0'}});
riskContent.add(loadingLabel);

top10.evaluate(function(fc) {
  riskContent.remove(loadingLabel);
  if (!fc || !fc.features || fc.features.length === 0) {
    riskContent.add(ui.Label({value:'No data — check connection',style:{fontSize:'11px',color:'#888888'}}));
    return;
  }
  fc.features.forEach(function(f, i) {
    var name = f.properties.name;
    var score = f.properties.risk_score || 0;
    var pct = (score * 100).toFixed(1);
    var barWidth = Math.round(score * 120);
    var dotColor = '#FF00FF';
    var textColor = i === 0 ? '#8B0000' : i < 3 ? '#FF0000' : '#333333';
    var displayName;
    var isUnknown = !name;
    if (isUnknown && f.geometry && f.geometry.coordinates && f.geometry.coordinates[0]) {
      var ring = f.geometry.coordinates[0];
      var sumLng = 0, sumLat = 0;
      for (var k = 0; k < ring.length; k++) { sumLng += ring[k][0]; sumLat += ring[k][1]; }
      var cLng = (sumLng / ring.length).toFixed(4);
      var cLat = (sumLat / ring.length).toFixed(4);
      displayName = cLat + '°N, ' + Math.abs(parseFloat(cLng)).toFixed(4) + '°W ▲ unregistered';
    } else {
      displayName = name || 'Unknown';
    }
    var deployAccess = f.properties.deploy_access || 0;
    var deployIcon = deployAccess > 0.5 ? '✓' : '–';
    var deployColor = deployAccess > 0.5 ? '#1B5E20' : '#BBBBBB';
    var basinName = f.properties.basin_name || '';
    var basinLabel = basinName ? ' · ' + basinName : '';
    var row = ui.Panel({
      widgets:[
        ui.Label({value:(i+1)+'.',style:{fontSize:'11px',color:textColor,fontWeight:'bold',margin:'0 4px 0 0',width:'18px'}}),
        ui.Label({value:displayName + basinLabel,style:{fontSize:'11px',color: isUnknown ? '#8B0000' : textColor,margin:'0',stretch:'horizontal',fontWeight:i===0?'bold':'normal',fontStyle: isUnknown ? 'italic' : 'normal'}}),
        ui.Label({value:pct+'%',style:{fontSize:'11px',color:dotColor,fontWeight:'bold',margin:'0 0 0 4px'}}),
        ui.Label({value:deployIcon,style:{fontSize:'13px',color:deployColor,fontWeight:'bold',margin:'0 0 0 6px'}})
      ],
      layout:ui.Panel.Layout.Flow('horizontal'),
      style:{margin:'0 0 1px 0',stretch:'horizontal'}
    });
    var barBg = ui.Panel({style:{backgroundColor:'#f0f0f0',margin:'0 0 4px 22px',height:'4px',stretch:'horizontal'}});
    var barFill = ui.Panel({style:{backgroundColor:dotColor,margin:'0',height:'4px',width:barWidth+'px'}});
    barBg.add(barFill);
    riskContent.add(row);
    riskContent.add(barBg);
  });
  riskContent.add(ui.Label({value:'▲ Unregistered — GPS coordinates for field validation. GalaSat-identified only.',style:{fontSize:'9px',color:'#8B0000',margin:'6px 0 0 0',fontStyle:'italic'}}));
  riskContent.add(ui.Label({value:'✓ = remediation equipment accessible | – = terrain restricts access',style:{fontSize:'9px',color:'#555555',margin:'2px 0 0 0'}}));
});

var riskToggleBtn = ui.Button({
  label:'► Biomonitoring Priority — Top 10 Communities at Risk',
  style:{fontWeight:'bold',fontSize:'11px',color:'#8B0000',backgroundColor:'white',border:'none',margin:'0 0 4px 0',stretch:'horizontal'},
  onClick:function() {
    riskPanelVisible = !riskPanelVisible;
    riskContent.style().set('shown', riskPanelVisible);
    riskToggleBtn.setLabel(riskPanelVisible ? '▼ Biomonitoring Priority — Top 10 Communities at Risk' : '► Biomonitoring Priority — Top 10 Communities at Risk');
  }
});
var riskPanel = ui.Panel({widgets:[riskToggleBtn, riskContent],style:{position:'bottom-left',padding:'8px 14px',backgroundColor:'white',width:'360px'}});
Map.add(riskPanel);

// ============================================================
// COMMUNITY SEARCH PANEL
// ============================================================
var searchPanel = ui.Panel({
  style:{position:'top-left',padding:'8px 14px',backgroundColor:'white',width:'260px'}
});
searchPanel.add(ui.Label({value:'Search Communities',style:{fontWeight:'bold',fontSize:'11px',color:'#8B0000',margin:'0 0 4px 0'}}));

var searchBox = ui.Textbox({placeholder:'Loading communities...',style:{stretch:'horizontal',margin:'0 0 0 0',fontSize:'11px'}});
var searchBtn = ui.Button({label:'Search',style:{fontSize:'11px',color:'#8B0000',margin:'0 0 2px 4px'}});
var searchRow = ui.Panel({widgets:[searchBox, searchBtn],layout:ui.Panel.Layout.Flow('horizontal'),style:{margin:'0 0 2px 0',stretch:'horizontal'}});
var hintLabel = ui.Label({value:'Press Enter or click Search',style:{fontSize:'9px',color:'#888888',margin:'0 0 4px 0'}});
var dropdownPanel = ui.Panel({style:{margin:'0',padding:'0',shown:false,backgroundColor:'white',border:'1px solid #cccccc'}});
var statusLabel = ui.Label({value:'Loading communities — please wait...',style:{fontSize:'9px',color:'#FF6D00',margin:'2px 0 0 0'}});
var communityCache = null;

function flyToLocation(lng, lat, name) {
  Map.setCenter(lng, lat, 14);
  searchBox.setValue(name);
  dropdownPanel.style().set('shown', false);
  statusLabel.setValue('Zoomed to: ' + name);
  statusLabel.style().set('color','#1B5E20');
}

function buildDropdown(text) {
  dropdownPanel.clear();
  if (text.length < 1) { dropdownPanel.style().set('shown', false); return; }
  if (!communityCache) {
    statusLabel.setValue('Still loading — try again in a moment...');
    statusLabel.style().set('color','#FF6D00');
    dropdownPanel.style().set('shown', false);
    return;
  }
  var lowerText = text.toLowerCase();
  var matches = [];
  for (var i = 0; i < communityCache.length; i++) {
    var c = communityCache[i];
    if (c.name && c.name.toLowerCase().indexOf(lowerText) === 0) matches.push(c);
  }
  if (matches.length === 0) {
    statusLabel.setValue('No matches for "' + text + '"');
    statusLabel.style().set('color','#888888');
    dropdownPanel.style().set('shown', false);
    return;
  }
  var showCount = Math.min(matches.length, 8);
  statusLabel.setValue(matches.length + ' match' + (matches.length > 1 ? 'es' : '') + ' — click to zoom');
  statusLabel.style().set('color','#1B5E20');
  for (var j = 0; j < showCount; j++) {
    (function(community) {
      var btn = ui.Button({
        label: community.name + (community.basin ? ' · ' + community.basin : ''),
        style:{stretch:'horizontal',fontSize:'11px',color:'#333333',backgroundColor:'white',border:'none',margin:'0',padding:'4px 8px',textAlign:'left'},
        onClick: function() { flyToLocation(community.lng, community.lat, community.name); }
      });
      dropdownPanel.add(btn);
    })(matches[j]);
  }
  if (matches.length > 8) dropdownPanel.add(ui.Label({value:'+ ' + (matches.length - 8) + ' more — keep typing to filter',style:{fontSize:'9px',color:'#888888',margin:'2px 4px 2px 8px'}}));
  dropdownPanel.style().set('shown', true);
}

exposedCommunities.evaluate(function(fc) {
  if (!fc || !fc.features) return;
  communityCache = [];
  fc.features.forEach(function(f) {
    var name = f.properties.name;
    var basin = f.properties.basin_name || '';
    if (!name) return;
    var coords = f.geometry ? f.geometry.coordinates : null;
    if (!coords) return;
    var lng, lat;
    if (typeof coords[0] === 'number') { lng = coords[0]; lat = coords[1]; }
    else if (coords[0] && typeof coords[0][0] === 'number') { lng = coords[0][0]; lat = coords[0][1]; }
    else if (coords[0] && coords[0][0]) {
      var ring = coords[0]; var sLng = 0, sLat = 0;
      for (var k = 0; k < ring.length; k++) { sLng += ring[k][0]; sLat += ring[k][1]; }
      lng = sLng / ring.length; lat = sLat / ring.length;
    }
    if (lng && lat) communityCache.push({name:name, basin:basin, lng:lng, lat:lat});
  });
  communityCache.sort(function(a,b){ return a.name.localeCompare(b.name); });
  statusLabel.setValue(communityCache.length + ' communities ready — type + press Enter or Search');
  statusLabel.style().set('color','#1B5E20');
  searchBox.setPlaceholder('Type name — e.g. Akomfere, Tano...');

});

function runSearch() {
  var text = searchBox.getValue();
  dropdownPanel.clear();
  dropdownPanel.style().set('shown', false);
  if (!text || text.length === 0) {
    statusLabel.setValue(communityCache ? communityCache.length + ' communities ready' : 'Loading...');
    statusLabel.style().set('color', communityCache ? '#1B5E20' : '#FF6D00');
    return;
  }
  if (!communityCache) { statusLabel.setValue('Still loading — try again...'); statusLabel.style().set('color','#FF6D00'); return; }
  statusLabel.setValue('Searching...'); statusLabel.style().set('color','#FF6D00');
  buildDropdown(text);
}

searchBox.onChange(function(text) { runSearch(); });
searchBtn.onClick(function() { runSearch(); });
searchPanel.add(searchRow);
searchPanel.add(hintLabel);
searchPanel.add(statusLabel);
searchPanel.add(dropdownPanel);

// ============================================================
// ALERT PANEL
// ============================================================
var alertVisible = false;
var alertContent = ui.Panel({style:{margin:'0',padding:'0',shown:false}});
var alertStatus = ui.Label({value:"Checking satellite data...",style:{fontSize:"11px",color:"#555555",margin:"0 0 4px 0"}});
alertContent.add(alertStatus);
newMiningHa.evaluate(function(ha) {
  if (ha && ha >= 0.1) {
    alertStatus.setValue("New mining detected this week");
    alertStatus.style().set("color","#8B0000");
    alertStatus.style().set("fontWeight","bold");
    alertContent.add(ui.Label({value:"Area: " + ha + " hectares",style:{fontSize:"14px",fontWeight:"bold",color:"#FF0000",margin:"2px 0"}}));
    alertContent.add(ui.Label({value:"Isolated sites: 0.5 ha min | Connected expansion: 0.1 ha min",style:{fontSize:"10px",color:"#555555",margin:"0 0 2px 0"}}));
    alertContent.add(ui.Label({value:"Toggle: Mining — New Activity This Week (red layer)",style:{fontSize:"10px",color:"#555555",margin:"0"}}));
  } else {
    alertStatus.setValue("No new mining activity detected this week.");
    alertStatus.style().set("color","#1B5E20");
    alertStatus.style().set("fontWeight","bold");
    alertContent.add(ui.Label({value:"Isolated: 0.5 ha min | Connected expansion: 0.1 ha min",style:{fontSize:"10px",color:"#888888",margin:"2px 0 0 0"}}));
    alertContent.add(ui.Label({value:"Recheck: next Sentinel-2 pass",style:{fontSize:"10px",color:"#888888",margin:"0"}}));
  }
});
var alertToggleBtn = ui.Button({
  label:'► ⚠ GalaSat Alert — New Activity',
  style:{fontWeight:'bold',fontSize:'12px',color:'#8B0000',backgroundColor:'white',border:'none',margin:'0 0 4px 0',stretch:'horizontal'},
  onClick:function(){
    alertVisible = !alertVisible;
    alertContent.style().set('shown', alertVisible);
    alertToggleBtn.setLabel(alertVisible ? '▼ ⚠ GalaSat Alert — New Activity' : '► ⚠ GalaSat Alert — New Activity');
  }
});
var alertPanel = ui.Panel({widgets:[alertToggleBtn, alertContent],style:{position:'top-center',padding:'8px 14px',backgroundColor:'white',width:'320px'}});
Map.add(alertPanel);

// ============================================================
// STATS PANEL
// ============================================================
var statsVisible = false;
var statsContent = ui.Panel({style:{margin:"0",padding:"0",shown:false}});
var statsPanel_outer = ui.Panel({style:{position:"top-right",padding:"8px 14px",backgroundColor:"white",width:"360px"}});
var statsToggleBtn = ui.Button({
  label:"► GalaSat v2.1 — Key Findings",
  style:{fontWeight:"bold",fontSize:"12px",color:"#8B0000",backgroundColor:"white",border:"none",margin:"0 0 4px 0",stretch:"horizontal"},
  onClick:function(){
    statsVisible = !statsVisible;
    statsContent.style().set("shown", statsVisible);
    statsToggleBtn.setLabel(statsVisible ? "▼ GalaSat v2.1 — Key Findings" : "► GalaSat v2.1 — Key Findings");
  }
});
statsPanel_outer.add(statsToggleBtn);
statsContent.add(ui.Label({value:"GalaSat v2.1 — Key Findings",style:{fontWeight:"bold",fontSize:"12px",color:"#8B0000",margin:"0 0 6px 0"}}));
var makeStatRow = function(label, value, color) {
  return ui.Panel({
    widgets:[
      ui.Label({value:label,style:{fontSize:"11px",color:"#333333",margin:"0 8px 0 0",stretch:"horizontal"}}),
      ui.Label({value:value,style:{fontSize:"11px",color:color||"#FF6D00",fontWeight:"bold",margin:"0"}})
    ],
    layout:ui.Panel.Layout.Flow("horizontal"),
    style:{margin:"0 0 4px 0",stretch:"horizontal"}
  });
};
statsContent.add(makeStatRow("Cumulative disturbance 2014-2026:","170,960 ha","#FF6D00"));
statsContent.add(makeStatRow("Direct spread zone (500m buffer):","821,711 ha","#D500F9"));
statsContent.add(makeStatRow("Ghana official estimate:","5,500 ha","#888888"));
statsContent.add(makeStatRow("Turbid water (wet season):","26,026 ha","#8B0000"));
statsContent.add(makeStatRow("Contaminated river segments:","52 confirmed (v2.1 AOI)","#8B0000"));
statsContent.add(makeStatRow("Communities within 2km — GalaSat (GMW+training, conservative floor):","310","#F50057"));
statsContent.add(makeStatRow("Communities within 2km — GMW:","208","#888888"));
statsContent.add(makeStatRow("Population within 2km of mining:","2,803,545","#F50057"));
statsContent.add(makeStatRow("Classifier accuracy:","93–95% | Kappa 0.87–0.91","#FF6D00"));
statsContent.add(ui.Label({value:"Safe deployment window: November — February",style:{fontSize:"10px",color:"#1B5E20",margin:"4px 0 0 0",fontWeight:"bold"}}));
statsPanel_outer.add(statsContent);
Map.add(statsPanel_outer);

// ============================================================
// LEGEND
// ============================================================
var legendVisible = false;
var legendContent = ui.Panel({style:{margin:'0',padding:'0',shown:false}});
var makeRow = function(color, label){
  return ui.Panel({widgets:[
    ui.Label({style:{backgroundColor:color,padding:'10px',margin:'0 6px 3px 0',border:'1px solid #cccccc'}}),
    ui.Label({value:label,style:{margin:'0 0 3px 6px',fontSize:'11px'}})
  ],layout:ui.Panel.Layout.Flow('horizontal')});
};
var makeSection = function(text){
  return ui.Label({value:text,style:{fontWeight:'bold',fontSize:'10px',margin:'5px 0 2px 0',color:'#8B0000'}});
};
legendContent.add(ui.Label({value:'Ashanti Region, Ghana | Remediation Sequencing + Biomonitoring Support',style:{fontSize:'10px',margin:'0 0 5px 0',color:'#555555'}}));
legendContent.add(makeSection('— MINING DISTURBANCE —'));
legendContent.add(makeRow('#FF6D00','Mining Activity — all years'));
legendContent.add(makeRow('#FFFF00','Mining Footprints — Global Mining Watch'));
legendContent.add(makeRow('#7B1FA2','AI-Detected Mining Disturbance 2025'));
legendContent.add(makeRow('#b71c1c','Iron Oxide Anomaly — ASGM Soil Contamination Signal'));
legendContent.add(makeRow('#FF0000','New Mining Activity This Week'));
legendContent.add(makeSection('— WATER CONTAMINATION —'));
legendContent.add(makeRow('#8B0000','Mercury-Contaminated Water (turbid, wet season)'));
legendContent.add(makeRow('#651FFF','Waterways Adjacent to Mining'));
legendContent.add(makeRow('#FFFFFF','Uncontaminated Water Bodies'));
legendContent.add(makeRow('#00BCD4','Seasonal Flood Zone — Deployment Safety'));
legendContent.add(makeSection('— MERCURY TRANSPORT —'));
legendContent.add(makeRow('#D500F9','Direct Spread Zone — 500m from Mining'));
legendContent.add(makeRow('#00897B','Downstream Low-Elevation Zones'));
legendContent.add(makeRow('#F9A825','Mercury Transport Pathways (SRTM-derived)'));
legendContent.add(makeRow('#7B1FA2','River Transport Risk — Downstream Network'));
legendContent.add(makeSection('— COMMUNITY EXPOSURE —'));
legendContent.add(makeRow('#01579b','Mercury Exposure Risk Score'));
legendContent.add(makeRow('#FF00FF','Biomonitoring Critical Risk — Ranks 1-10'));
legendContent.add(makeRow('#00E5FF','Biomonitoring High Risk — Ranks 11-50'));
legendContent.add(makeRow('#0D47A1','Biomonitoring Elevated Risk — Ranks 51-100'));
legendContent.add(makeRow('#76FF03','All Communities within 2km of Mining'));
legendContent.add(makeRow('#FF6F00','Population within 2km of Mining (WorldPop)'));
legendContent.add(makeRow('#FFF9C4','Night Lights 2024 — Economic Activity (VIIRS)'));
legendContent.add(makeSection('— DEPLOYMENT & ACCESS —'));
legendContent.add(makeRow('#FF3D00','Borehole Exclusion Zone — Groundwater Contamination Risk'));
legendContent.add(makeRow('#00BCD4','Seasonal Flood Zone — Deployment Safety'));
legendContent.add(makeRow('#76FF03','Remediation-Accessible Terrain'));
legendContent.add(makeSection('— FARMING SAFETY —'));
legendContent.add(makeRow('#FF1493','Cropland at Mercury Risk — Food Chain Contamination'));
legendContent.add(makeRow('#FFFFFF','Cropland Safe from Contamination'));
legendContent.add(makeSection('— ENVIRONMENT —'));
legendContent.add(makeRow('#1B5E20','Forest Cover'));
legendContent.add(makeRow('#4575b4','Wet Season Precipitation Jun-Sep (GPM)'));
legendContent.add(makeRow('#fed98e','Dry Season Precipitation Nov-Feb (GPM)'));
legendContent.add(makeSection('— REFERENCE —'));
legendContent.add(makeRow('#0091EA','HydroSHEDS Rivers'));
legendContent.add(makeRow('#FFD600','HydroSHEDS Watersheds'));
legendContent.add(makeRow('#616161','Major Roads'));
legendContent.add(makeRow('#B0BEC5','River Basin Boundaries'));
legendContent.add(makeRow('#90A4AE','Extended Cumulative Mining 1984-2026'));
legendContent.add(makeRow('#FFFFFF','Pilot AOI Boundary'));
legendContent.add(ui.Label({value:'↑ N',style:{fontWeight:'bold',fontSize:'12px',color:'#333333',margin:'5px 0 2px 0'}}));
legendContent.add(ui.Label({value:'93-95% accuracy | Kappa 0.87-0.91 | 628 hand-mapped polygons | pre-grant build',style:{fontSize:'9px',color:'#555555',margin:'4px 0 1px 0'}}));
legendContent.add(ui.Label({value:'Validated: Oct 2025–Feb 2026 | EPSG:4326 | SCL+QA60 cloud mask',style:{fontSize:'9px',color:'#555555',margin:'0 0 1px 0'}}));
legendContent.add(ui.Label({value:'Sources: GMW | JRC | USGS SRTM | HydroSHEDS | ESA WorldCover | Ghana Statistical Service 2021',style:{fontSize:'9px',color:'#555555',margin:'0 0 1px 0'}}));
legendContent.add(ui.Label({value:'Generated: June 2026 | INFN8 VZN | restoreghana.ca | info@infn8vzn.org',style:{fontSize:'9px',color:'#8B0000',margin:'0'}}));
var toggleBtn = ui.Button({
  label: '► GalaSat v2.1 — Mercury Contamination Intelligence',
  style:{fontWeight:'bold',fontSize:'12px',color:'#8B0000',backgroundColor:'white',border:'none',margin:'0 0 4px 0',stretch:'horizontal'},
  onClick: function() {
    legendVisible = !legendVisible;
    legendContent.style().set('shown', legendVisible);
    toggleBtn.setLabel(legendVisible ? '▼ GalaSat v2.1 — Mercury Contamination Intelligence' : '► GalaSat v2.1 — Mercury Contamination Intelligence');
  }
});
var legend = ui.Panel({widgets:[toggleBtn, legendContent],style:{position:'bottom-right',padding:'8px 14px',backgroundColor:'white',maxHeight:'540px'}});
Map.add(legend);


// TILE LOADING INDICATOR — bottom-center floating label
// Appears when layers are toggled, clears after 10 seconds
var tileLoadingLabel = ui.Label({
  value:'',
  style:{fontSize:'12px',fontWeight:'bold',color:'white',backgroundColor:'#8B0000',padding:'6px 16px',margin:'0'}
});
var tileLoadingPanel = ui.Panel({
  widgets:[tileLoadingLabel],
  style:{position:'bottom-center',padding:'0',backgroundColor:'rgba(0,0,0,0)',shown:false}
});
Map.add(tileLoadingPanel);

var loadingTimer = null;
function showTileLoading(layerName) {
  tileLoadingLabel.setValue('⟳ Loading tiles: ' + layerName + ' — please wait...');
  tileLoadingPanel.style().set('shown', true);
  // Clear after 10 seconds
  ui.util.setTimeout(function() {
    tileLoadingPanel.style().set('shown', false);
  }, 10000);
}

// layer.onChange not supported in GEE — loading indicator triggered by map clicks instead
Map.onClick(function() { showTileLoading('map'); });

// Search panel added last — renders on top
Map.add(searchPanel);

// ============================================================
// EXPORTS
// ============================================================
Export.image.toDrive({image:rfMining,description:'INFN8VZN_RF_Mining_2025',folder:'GalaSat',fileNamePrefix:'infn8vzn_rf_mining_2025',region:AOI,scale:20,crs:'EPSG:4326',maxPixels:1e13,fileFormat:'GeoTIFF'});
Export.image.toDrive({image:turbidWater,description:'INFN8VZN_Turbid_Water_2025',folder:'GalaSat',fileNamePrefix:'infn8vzn_turbid_water_2025',region:AOI,scale:20,crs:'EPSG:4326',maxPixels:1e13,fileFormat:'GeoTIFF'});
Export.image.toDrive({image:floodZone,description:'INFN8VZN_Flood_Zone',folder:'GalaSat',fileNamePrefix:'infn8vzn_flood_zone',region:AOI,scale:30,crs:'EPSG:4326',maxPixels:1e13,fileFormat:'GeoTIFF'});
Export.image.toDrive({image:cumulativeMining,description:'INFN8VZN_Cumulative_Mining_2014_2026',folder:'GalaSat',fileNamePrefix:'infn8vzn_cumulative_mining_2014_2026',region:AOI,scale:30,crs:'EPSG:4326',maxPixels:1e13,fileFormat:'GeoTIFF'});
Export.image.toDrive({image:mercuryFlowRisk,description:'INFN8VZN_Mercury_Transport_Pathways',folder:'GalaSat',fileNamePrefix:'infn8vzn_mercury_transport_pathways',region:AOI,scale:30,crs:'EPSG:4326',maxPixels:1e13,fileFormat:'GeoTIFF'});
Export.image.toDrive({image:deploymentAccess,description:'INFN8VZN_Remediation_Accessible_Terrain',folder:'GalaSat',fileNamePrefix:'infn8vzn_remediation_accessible_terrain',region:AOI,scale:30,crs:'EPSG:4326',maxPixels:1e13,fileFormat:'GeoTIFF'});
Export.image.toDrive({image:newMiningImg,description:'INFN8VZN_New_Mining_Activity_This_Week',folder:'GalaSat',fileNamePrefix:'infn8vzn_new_mining_this_week',region:AOI,scale:20,crs:'EPSG:4326',maxPixels:1e13,fileFormat:'GeoTIFF'});
Export.image.toDrive({image:exposedPop,description:'INFN8VZN_Population_Exposure_2km',folder:'GalaSat',fileNamePrefix:'infn8vzn_population_exposure_2km',region:AOI,scale:100,crs:'EPSG:4326',maxPixels:1e13,fileFormat:'GeoTIFF'});
Export.image.toDrive({image:gpmWetSeason,description:'INFN8VZN_GPM_Wet_Season_Precipitation',folder:'GalaSat',fileNamePrefix:'infn8vzn_gpm_wet_season',region:AOI,scale:1000,crs:'EPSG:4326',maxPixels:1e13,fileFormat:'GeoTIFF'});
Export.image.toDrive({image:directSpreadZone,description:'INFN8VZN_Direct_Spread_Zone_500m',folder:'GalaSat',fileNamePrefix:'infn8vzn_direct_spread_zone_500m',region:AOI,scale:30,crs:'EPSG:4326',maxPixels:1e13,fileFormat:'GeoTIFF'});
Export.image.toDrive({image:highRiverRisk,description:'INFN8VZN_River_Transport_Risk',folder:'GalaSat',fileNamePrefix:'infn8vzn_river_transport_risk',region:AOI,scale:30,crs:'EPSG:4326',maxPixels:1e13,fileFormat:'GeoTIFF'});
Export.image.toDrive({image:boreholeExclusionZone,description:'INFN8VZN_Borehole_Exclusion_Zone',folder:'GalaSat',fileNamePrefix:'infn8vzn_borehole_exclusion_zone',region:AOI,scale:1000,crs:'EPSG:4326',maxPixels:1e13,fileFormat:'GeoTIFF'});
Export.image.toDrive({image:cropAtRisk,description:'INFN8VZN_Cropland_At_Risk',folder:'GalaSat',fileNamePrefix:'infn8vzn_cropland_at_risk',region:AOI,scale:500,crs:'EPSG:4326',maxPixels:1e13,fileFormat:'GeoTIFF'});
Export.image.toDrive({image:cropSafe,description:'INFN8VZN_Cropland_Safe',folder:'GalaSat',fileNamePrefix:'infn8vzn_cropland_safe',region:AOI,scale:500,crs:'EPSG:4326',maxPixels:1e13,fileFormat:'GeoTIFF'});
Export.image.toDrive({image:worldCover.eq(40).selfMask().rename('cropland'),description:'INFN8VZN_WorldCover_Cropland',folder:'GalaSat',fileNamePrefix:'infn8vzn_worldcover_cropland',region:AOI,scale:100,crs:'EPSG:4326',maxPixels:1e13,fileFormat:'GeoTIFF'});

print('GalaSat v2.1 — INFN8 VZN | COMPLETE');
print('Mining detector: improved — 3 new conditions for bare excavated earth');
print('New activity: always-on bright red toggle + alert panel');
print('Legend: collapsible — click title to open/close');
print('Stats panel: top-right key figures on map');
print('SRTM: elevation, slope, mercury transport, deployment access');
print('Iron oxide anomaly: threshold 2.5 — ASGM soil contamination signal');
print('Community dots: 4-tier — Fuchsia top 10 | Cyan 11-50 | Deep Blue 51-100 | Lime all remaining 2km communities');
print('Farming Safety: WorldCover cropland crossed with contamination envelope');
print('Borehole Exclusion Zone: GLOBGM shallow water table + contamination proximity');
print('Layer status panel: 12 tracked computations with percentage complete');
print('Search: starts-with autocomplete, 310 communities, basin names');
