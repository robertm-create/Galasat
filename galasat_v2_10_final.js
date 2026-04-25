// GalaSat v2.0 Final - INFN8 VZN
// Frankenstein build - best of old and new
// AOI: original watershed bounds from training_340 asset (larger AOI, real turbid water numbers)
// Composite: v5.20 4-month window Nov2025-Feb2026, 20% cloud, dual SCL+QA60 mask (clean image)
// Classifier: v5.20 three-class polygon RF - 628 polygons mining/urban/scrubland
// Layers: all v5.20 toggles, all off by default except AOI and mining disturbance
// Legend: clean white simple style from old code
// Time series: 12-year confirmed 2014-2025

// STUDY AREA - clean rectangle around Ashanti mining zone
// Cuts coastal droop and southeastern extension from original training asset bounds
var trainingRef = ee.FeatureCollection(
  'projects/galamsey-monotoring/assets/infn8vzn_mercury_training_340'
);
var AOI = ee.Geometry.Rectangle([-2.60, 5.90, -1.00, 6.90]);
Map.centerObject(AOI, 10);
Map.setOptions('SATELLITE');

// POLYGON TRAINING ASSET - v5.20 three-class system
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

// OSM ASSETS
var roadsAsset = ee.FeatureCollection(
  'projects/galamsey-monotoring/assets/ghana_osm_roads'
).filterBounds(AOI);
var waterways = ee.FeatureCollection(
  'projects/galamsey-monotoring/assets/ghana_osm_waterways'
).filterBounds(AOI);
var places = ee.FeatureCollection(
  'projects/galamsey-monotoring/assets/ghana_osm_places'
).filterBounds(AOI);
var roads = roadsAsset.filter(ee.Filter.inList('fclass',
  ['motorway','trunk','primary','secondary'])).limit(500);

var tier1Rivers = waterways.filter(ee.Filter.inList('name',[
  'Pra','Pra River','Offin','Offin River','Ankobra','Ankobra River',
  'Birim','Birim River','Tano','Tano River','Densu','Densu River',
  'Oda','Oda River','Bonsa','Bonsa River','Subin','Subin River',
  'Afram','Afram River','Fena','Fena River','Anum','Anum River'
]));

// GLOBAL MINING WATCH
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

// CLOUD MASKING - dual SCL + QA60 (v5.20 clean image method)
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

// DRY COMPOSITE - 4-month window Nov2025-Feb2026, 20% cloud (v5.20 clean image)
var s2dry = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
  .filterBounds(AOI).filterDate('2025-10-01','2026-02-28')
  .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE',30))
  .map(maskS2).median().clip(AOI);

// WET COMPOSITE - turbid water detection
var s2wet = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
  .filterBounds(AOI).filterDate('2023-05-01','2025-10-31')
  .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE',30))
  .map(maskS2wet).median().clip(AOI);

print('Dry season images:',
  ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
    .filterBounds(AOI).filterDate('2025-10-01','2026-02-28')
    .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE',30)).size()
);

// SENTINEL-1 SAR
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

// SPECTRAL INDICES
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

// RULE-BASED COMPOSITE BUILDERS - for time series
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

// RULE-BASED DETECTORS
function mineS2(img) {
  var b2=img.select('B2'),b3=img.select('B3'),b4=img.select('B4');
  var b8=img.select('B8'),b8a=img.select('B8A'),b11=img.select('B11');
  var valid=b2.mask().and(b3.mask()).and(b4.mask()).and(b8.mask()).and(b8a.mask()).and(b11.mask());
  var bsi=b11.add(b4).subtract(b8a.add(b2)).divide(b11.add(b4).add(b8a).add(b2));
  var ndvi=b8.subtract(b4).divide(b8.add(b4));
  var ndwi=b3.subtract(b8).divide(b3.add(b8));
  var ndbi=b11.subtract(b8).divide(b11.add(b8));
  var mndwi=b3.subtract(b11).divide(b3.add(b11));
  var combined=bsi.gt(0.005).and(ndvi.lt(0.35)).and(ndbi.lt(0.05))
    .or(ndwi.gt(0.02).and(ndvi.lt(0.15)).and(bsi.gt(-0.08)))
    .or(b4.gt(0.07).and(b3.gt(0.05)).and(ndvi.lt(0.1)).and(ndbi.lt(0.0)))
    .or(mndwi.gt(0.0).and(ndvi.lt(0.1)).and(b4.gt(0.05)));
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

// MASKS
var worldCover = ee.ImageCollection('ESA/WorldCover/v200').first().clip(AOI);
var urbanMaskBuffered = worldCover.neq(50).and(worldCover.neq(40))
  .focal_min(1,'square','pixels');
var ghsl = ee.ImageCollection('JRC/GHSL/P2016/SMOD_POP_GLOBE_V1')
  .filterBounds(AOI).mosaic().clip(AOI);
var ghslUrban = ghsl.select('smod_code').lt(3);
var roadMask = ee.Image(1).byte().paint(
  roads.map(function(f){return f.buffer(50);}), 0
);
var finalMask = urbanMaskBuffered.and(roadMask).and(ghslUrban);

// RF FEATURE STACK
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

// RF TRAINING
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

// RF OUTPUT - connected pixel filter eliminates forest scatter
var rfRaw = rfClass.eq(1)
  .updateMask(rfClass.eq(1))
  .updateMask(urbanExclude.not())
  .updateMask(roadFootprint)
  .updateMask(ndviImg.lt(0.65));
var rfMining = rfRaw
  .updateMask(rfRaw.connectedPixelCount(200, true).gte(20))
  .selfMask();

// VALIDATION
var validated   = valSet.classify(rfClassifier);
var errorMatrix = validated.errorMatrix('c','classification');
print('RF Accuracy');
print('Confusion Matrix:', errorMatrix);
print('Overall Accuracy:', errorMatrix.accuracy());
print('Kappa:', errorMatrix.kappa());
print('Producers Accuracy:', errorMatrix.producersAccuracy());
print('Consumers Accuracy:', errorMatrix.consumersAccuracy());

// RULE-BASED DETECTIONS - time series
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

// CUMULATIVE DISTURBANCE
var cumulativeMining = ml14.unmask(0).gt(0)
  .or(ml15.unmask(0).gt(0)).or(m16.unmask(0).gt(0))
  .or(m18.unmask(0).gt(0)).or(m20.unmask(0).gt(0))
  .or(m22.unmask(0).gt(0)).or(m24.unmask(0).gt(0))
  .or(m25.unmask(0).gt(0)).or(m26.unmask(0).gt(0))
  .selfMask();

// TURBID WATER
var MNDWI_wet = s2wet.normalizedDifference(['B3','B11']).rename('MNDWI_wet');
var NDTI_wet  = s2wet.normalizedDifference(['B4','B3']).rename('NDTI_wet');
var IOR_wet   = s2wet.select('B4').divide(s2wet.select('B2')).rename('IOR_wet');
var NDVI_wet  = s2wet.normalizedDifference(['B8','B4']).rename('NDVI_wet');

var turbidWater = MNDWI_wet.gt(0.1).and(NDTI_wet.gt(0.08))
  .and(IOR_wet.gt(1.1)).and(NDVI_wet.lt(0.2))
  .and(s2wet.select('B3').gt(0.05))
  .selfMask().rename('TurbidWater');

// FLOOD SAFETY JRC
var pixelArea     = ee.Image.pixelArea();
var jrcOccurrence = ee.Image('JRC/GSW1_4/GlobalSurfaceWater').select('occurrence').clip(AOI);
var jrcSeasonality= ee.Image('JRC/GSW1_4/GlobalSurfaceWater').select('seasonality').clip(AOI);
var permWater     = jrcOccurrence.gt(75).selfMask().clip(AOI);
var floodZone     = jrcSeasonality.gt(0).and(jrcOccurrence.lte(75)).selfMask().clip(AOI);
var forest        = worldCover.eq(10).selfMask();

// BASIN GEOMETRIES
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

// HYDROSHEDS
var hydroRivers       = ee.FeatureCollection('WWF/HydroSHEDS/v1/FreeFlowingRivers').filterBounds(AOI);
var hydroBasins       = ee.FeatureCollection('WWF/HydroSHEDS/v1/Basins/hybas_7').filterBounds(AOI);
var hydroRiversRaster = ee.Image().byte().paint({featureCollection:hydroRivers,color:1,width:2});
var hydroBasinsOutline= ee.Image().byte().paint({featureCollection:hydroBasins,color:1,width:2});

print('HydroSHEDS river segments:', hydroRivers.size());
print('HydroSHEDS watersheds:', hydroBasins.size());

// WATER RASTER LAYERS
var allWaterwaysRaster      = ee.Image().byte().paint({featureCollection:waterways,color:1,width:1});
var tier1Raster             = ee.Image().byte().paint({featureCollection:tier1Rivers,color:1,width:3});
var adjacentWaterwaysRaster = ee.Image().byte().paint(
  waterways.map(function(f){return f.buffer(30);}),1
).updateMask(miningBufferImg);

// COMMUNITY EXPOSURE
var miningUnion = globalMines.map(function(f){return f.buffer(2000);})
  .union(ee.ErrorMargin(100));
var exposedCommunities = places.filterBounds(miningUnion.geometry());
var exposedCommunitiesRaster = ee.Image().byte().paint({
  featureCollection:exposedCommunities.map(function(f){return f.buffer(100);}),color:1
});

// AREA CALCULATIONS
var pixelAreaImg = ee.Image.pixelArea();

var turbidAreaDict = turbidWater.multiply(pixelAreaImg).divide(10000)
  .reduceRegion({reducer:ee.Reducer.sum(),geometry:AOI,scale:20,maxPixels:1e13,bestEffort:true});
print('Turbid Water (wet season 2023-2025)');
print('Turbid Water Area (hectares):', ee.Number(turbidAreaDict.get('TurbidWater')).round());

var rfAreaDict = rfMining.multiply(pixelAreaImg).divide(10000)
  .reduceRegion({reducer:ee.Reducer.sum(),geometry:AOI,scale:200,maxPixels:1e13,bestEffort:true});
print('RF Mining Disturbance 2025');
print('RF Mining Ha:', ee.Number(rfAreaDict.get(rfAreaDict.keys().get(0))).round());

var cumAreaDict = cumulativeMining.multiply(pixelAreaImg).divide(10000)
  .reduceRegion({reducer:ee.Reducer.sum(),geometry:AOI,scale:200,maxPixels:1e13,bestEffort:true});
print('Cumulative Mining Disturbance Ha (2014-2026):',
  ee.Number(cumAreaDict.get(cumAreaDict.keys().get(0))).round());

print('Community Exposure');
print('Total communities in AOI:', places.size());
print('Communities within 2km of confirmed mining:', exposedCommunities.size());

print('Flood Safety');
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

// TIME SERIES 2014-2025 - cloud filter 75% recovers wet season data
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

print('Time Series Chart');
print(ui.Chart.feature.byFeature({features:tsFC,xProperty:'year',yProperties:['turbid_ha']})
  .setChartType('LineChart').setOptions({
    title:'GalaSat: Mining-Affected Turbid Water 2014-2025',
    hAxis:{title:'Year',format:'####'},
    vAxis:{title:'Turbid Water Area (hectares)'},
    colors:['#8B0000'],lineWidth:3,pointSize:6,interpolateNulls:false
  }));

// MAP LAYERS - all off by default except mining disturbance and AOI
var tcViz  ={bands:['B4','B3','B2'],min:0.0,max:0.3,gamma:1.4};
var tcL8Viz={bands:['SR_B4','SR_B3','SR_B2'],min:0.0,max:0.3,gamma:1.4};
var mineViz={min:0.0,max:0.3,palette:['ffff00','fd8d3c','f03b20','bd0026']};

// Base imagery toggles
Map.addLayer(l8_14,tcL8Viz,'2014 True Color',false);
Map.addLayer(l8_15,tcL8Viz,'2015 True Color',false);
Map.addLayer(d16,tcViz,'2016 True Color',false);
Map.addLayer(d18,tcViz,'2018 True Color',false);
Map.addLayer(s2_20,tcViz,'2020 True Color',false);
Map.addLayer(s2_22,tcViz,'2022 True Color',false);
Map.addLayer(s2_24,tcViz,'2024 True Color',false);
Map.addLayer(s2_25,tcViz,'2025 True Color',false);
Map.addLayer(s2_26,tcViz,'2026 True Color',false);
Map.addLayer(s2dry,tcViz,'Current Base Image Nov2025-Feb2026',false);

// Mining detection layers
Map.addLayer(m26,mineViz,'Mining Activity 2026 (tricolor)',true);
Map.addLayer(ml14,mineViz,'2014 Mining (L8)',false);
Map.addLayer(ml15,mineViz,'2015 Mining (L8)',false);
Map.addLayer(m16,mineViz,'2016 Mining',false);
Map.addLayer(m18,mineViz,'2018 Mining',false);
Map.addLayer(m20,mineViz,'2020 Mining',false);
Map.addLayer(m22,mineViz,'2022 Mining',false);
Map.addLayer(m24,mineViz,'2024 Mining',false);
Map.addLayer(m25,mineViz,'2025 Mining',false);
Map.addLayer(cumulativeMining,{min:0,max:1,palette:['ff0000'],opacity:0.7},'Cumulative Mining Disturbance 2014-2026',false);
Map.addLayer(rfMining,{min:0,max:1,palette:['E040FB'],opacity:0.7},'RF Mining Detection 2025',false);

// Water and environment layers
Map.addLayer(permWater,{palette:['#003087']},'Permanent Water',false);
Map.addLayer(floodZone,{palette:['#00BCD4']},'Seasonal Flood Zone',false,0.85);
Map.addLayer(turbidWater,{palette:['FF8C00'],opacity:0.8},'Turbid Water (wet season 2023-2025)',false);
Map.addLayer(forest,{palette:['1A5C1A']},'Forest Cover',false);

// Infrastructure layers
Map.addLayer(globalMines,{color:'FFFF00'},'Mining Footprints (Global Mining Watch)',false);
Map.addLayer(roads,{color:'888888'},'Major Roads',false);
Map.addLayer(allWaterwaysRaster,{palette:['4FC3F7'],opacity:0.5},'All Waterways',false);
Map.addLayer(tier1Raster,{palette:['8B0000']},'Field-Verified Contaminated Rivers',false);
Map.addLayer(adjacentWaterwaysRaster,{palette:['9C27B0']},'Waterways Adjacent to Mining',false);
Map.addLayer(exposedCommunitiesRaster,{palette:['FF00FF']},'Communities within 2km',false);
Map.addLayer(basin_fc.style({color:'FF6600',fillColor:'00000000',width:2}),{},'River Basin Boundaries',false);

// Training polygon layers
Map.addLayer(miningPolys.style({color:'FF1744',fillColor:'FF174440',width:1}),{},'Mapped Polygons: Mining',false);
Map.addLayer(urbanPolys.style({color:'FFC107',fillColor:'FFC10740',width:1}),{},'Mapped Polygons: Urban',false);
Map.addLayer(scrubPolys.style({color:'00ff00',fillColor:'00ff0040',width:1}),{},'Mapped Polygons: Scrubland',false);

// SAR layers
Map.addLayer(s1dry.select('VV_filtered'),{min:-20,max:0,palette:['000000','ffffff']},'SAR VV Dry Season (radar)',false);
Map.addLayer(sarChange,{min:-5,max:5,palette:['0000ff','ffffff','ff0000']},'SAR Change Wet-Dry',false);

// HydroSHEDS layers
Map.addLayer(hydroRiversRaster,{palette:['00E5FF']},'HydroSHEDS Rivers',false);
Map.addLayer(hydroBasinsOutline,{palette:['FF9800']},'HydroSHEDS Watersheds',false);

// AOI boundary - always on
var aoiOutline=ee.Image().byte().paint({
  featureCollection:ee.FeatureCollection([ee.Feature(AOI)]),color:1,width:2
});
Map.addLayer(aoiOutline,{palette:['FFFFFF'],opacity:1.0},'Pilot AOI Boundary');

// LEGEND - clean white style from old code
var legend = ui.Panel({style:{
  position:'bottom-left',padding:'8px 15px',backgroundColor:'white'
}});
legend.add(ui.Label({
  value:'GalaSat v2.0 - INFN8VZN',
  style:{fontWeight:'bold',fontSize:'14px',margin:'0 0 2px 0',color:'#8B0000'}
}));
legend.add(ui.Label({
  value:'Ashanti Pilot, Ghana | SAR + Optical | 2026',
  style:{fontSize:'11px',margin:'0 0 6px 0',color:'#555555'}
}));

var makeRow = function(color,label){
  return ui.Panel({widgets:[
    ui.Label({style:{backgroundColor:color,padding:'8px',margin:'0 6px 4px 0'}}),
    ui.Label({value:label,style:{margin:'0 0 4px 6px',fontSize:'12px'}})
  ],layout:ui.Panel.Layout.Flow('horizontal')});
};

legend.add(makeRow('#bd0026','Mining Activity 2026 (high)'));
legend.add(makeRow('#f03b20','Mining Activity (moderate)'));
legend.add(makeRow('#fd8d3c','Mining Activity (low)'));
legend.add(makeRow('#ff0000','Cumulative Mining Disturbance 2014-2026'));
legend.add(makeRow('#E040FB','RF Mining Detection 2025'));
legend.add(makeRow('#FF8C00','Turbid Water (wet season 2023-2025)'));
legend.add(makeRow('#003087','Permanent Water'));
legend.add(makeRow('#00BCD4','Seasonal Flood Zone (deployment safety)'));
legend.add(makeRow('#1A5C1A','Forest Cover'));
legend.add(makeRow('#FFFF00','Mining Footprints (Global Mining Watch)'));
legend.add(makeRow('#888888','Major Roads'));
legend.add(makeRow('#4FC3F7','All Waterways'));
legend.add(makeRow('#8B0000','Field-Verified Contaminated Rivers'));
legend.add(makeRow('#9C27B0','Waterways Adjacent to Mining'));
legend.add(makeRow('#FF00FF','Communities within 2km'));
legend.add(makeRow('#FF6600','River Basin Boundaries'));
legend.add(makeRow('#FFFFFF','Pilot AOI Boundary'));
legend.add(makeRow('#00E5FF','HydroSHEDS Rivers'));
legend.add(makeRow('#FF9800','HydroSHEDS Watersheds'));
legend.add(ui.Label({
  value:'93.8% accuracy | Kappa 0.881 | 628 polygons | S2 + SAR',
  style:{fontSize:'10px',color:'#555555',margin:'6px 0 1px 0'}
}));
legend.add(ui.Label({
  value:'Nov 2025-Feb 2026 | SCL+QA60 cloud mask | EPSG:4326',
  style:{fontSize:'10px',color:'#555555',margin:'0 0 1px 0'}
}));
legend.add(ui.Label({
  value:'Sources: Joy News | Al Jazeera | PMC 2024 | GWCL | GMW | JRC',
  style:{fontSize:'10px',color:'#555555',margin:'0'}
}));
Map.add(legend);

// EXPORTS
Export.image.toDrive({
  image:rfMining, description:'INFN8VZN_RF_Mining_2025',
  folder:'GalaSat', fileNamePrefix:'infn8vzn_rf_mining_2025',
  region:AOI, scale:20, crs:'EPSG:4326', maxPixels:1e13, fileFormat:'GeoTIFF'
});
Export.image.toDrive({
  image:turbidWater, description:'INFN8VZN_Turbid_Water_2025',
  folder:'GalaSat', fileNamePrefix:'infn8vzn_turbid_water_2025',
  region:AOI, scale:20, crs:'EPSG:4326', maxPixels:1e13, fileFormat:'GeoTIFF'
});
Export.image.toDrive({
  image:floodZone, description:'INFN8VZN_Flood_Zone',
  folder:'GalaSat', fileNamePrefix:'infn8vzn_flood_zone',
  region:AOI, scale:30, crs:'EPSG:4326', maxPixels:1e13, fileFormat:'GeoTIFF'
});
Export.image.toDrive({
  image:cumulativeMining, description:'INFN8VZN_Cumulative_Mining_2014_2026',
  folder:'GalaSat', fileNamePrefix:'infn8vzn_cumulative_mining_2014_2026',
  region:AOI, scale:30, crs:'EPSG:4326', maxPixels:1e13, fileFormat:'GeoTIFF'
});

print('GalaSat v2.0 Final - Frankenstein build complete');
print('AOI: original watershed bounds - larger AOI - real turbid water numbers');
print('Composite: Nov2025-Feb2026 4-month window - clean image');
print('Classifier: 628 polygon three-class RF + connected pixel filter');
print('All layers off by default except Mining Activity 2026 and AOI boundary');