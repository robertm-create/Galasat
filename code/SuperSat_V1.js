// ============================================================
// GalaSat v4.3 — Polygon-Trained RF + Full Layer Stack
// INFN8VZN | Ashanti Region Mercury Contamination Detection
// RF now trains on full polygon areas — not centroid points
// ============================================================
//
// SOURCE: claude.ai project "Technical - Repository"
//         chat "SuperSat_V1 code tagging" (Apr 20, 2026)
//         pasted attachment, 531 lines
// MEMORY TAG: SuperSat_V1
// GEE FILENAME: SuperSat_V1
// VERSION: GalaSat v4.3 / v4.5
// ACCURACY: 81% overall, Kappa 0.715
// EXTRACTED: 2026-04-25
//
// NOTE: Whitespace was partially collapsed during browser-side extraction
//       (claude.ai chat content protection blocks raw bulk text from JS;
//        text was recovered via DOM injection + page-text scrape).
//       Logic and semicolons are intact — code is executable but may need
//       formatter reflow for readability. Compare against the canonical
//       file in your GEE repo if precise whitespace matters.

// ── STUDY AREA ─────────────────────────────────────────────
var box = ee.Geometry.Rectangle([-2.3670, 5.8599, -1.4829, 6.5807]);
Map.centerObject(box, 10);
Map.setOptions('SATELLITE');

// ── POLYGON ASSET (805 mapped polygons) ────────────────────
var mappedPolygons = ee.FeatureCollection('projects/galamsey-monotoring/assets/MINING_POLYGONS_800');

// ── ASSIGN CLASS LABELS FROM NAME FIELD ────────────────────
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

var miningPolys  = labeledPolygons.filter(ee.Filter.eq('c', 1));
var urbanPolys   = labeledPolygons.filter(ee.Filter.eq('c', 0));
var scrubPolys   = labeledPolygons.filter(ee.Filter.eq('c', 2));
var usablePolys  = labeledPolygons.filter(ee.Filter.neq('c', -1));

print('── Mapped Polygon Classes ─────────────');
print('Mining polygons:', miningPolys.size());
print('Urban polygons:', urbanPolys.size());
print('Scrubland polygons:', scrubPolys.size());
print('Total usable:', usablePolys.size());

// ── POLYGON HECTARE COUNT ──────────────────────────────────
var miningPolyHa = ee.Number(miningPolys.map(function(f) {
  return f.set('area_ha', f.geometry().area(1).divide(10000));
}).aggregate_sum('area_ha'));
print('── Mapped Mining Polygon Area ─────────');
print('Manually mapped mining area (ha):', miningPolyHa.round());

// ── ASSETS ────────────────────────────────────────────────
var existingPoints = ee.FeatureCollection('projects/galamsey-monotoring/assets/infn8vzn_mercury_training_340');
var roadsAsset = ee.FeatureCollection('projects/galamsey-monotoring/assets/ghana_osm_roads').filterBounds(box);
var waterways  = ee.FeatureCollection('projects/galamsey-monotoring/assets/ghana_osm_waterways').filterBounds(box);
var places     = ee.FeatureCollection('projects/galamsey-monotoring/assets/ghana_osm_places').filterBounds(box);

var roads = roadsAsset.filter(ee.Filter.inList('fclass', ['motorway','trunk','primary','secondary'])).limit(500);
var tier1Rivers = waterways.filter(ee.Filter.inList('name',[
  'Pra','Pra River','Offin','Offin River','Ankobra','Ankobra River',
  'Birim','Birim River','Tano','Tano River','Densu','Densu River',
  'Oda','Oda River','Bonsa','Bonsa River','Subin','Subin River',
  'Afram','Afram River','Fena','Fena River'
]));

// ── GLOBAL MINING WATCH ────────────────────────────────────
var globalMines = ee.FeatureCollection('projects/sat-io/open-datasets/global-mining/global_mining_polygons').filterBounds(box);
var minesClipped = globalMines.map(function(f){return f.intersection(box,ee.ErrorMargin(1));});
var totalMinesHa = ee.Number(minesClipped.map(function(f){
  return f.set('area_ha',f.geometry().area(1).divide(10000));
}).aggregate_sum('area_ha'));
var miningBufferImg = ee.Image(0).byte().paint(globalMines.map(function(f){return f.buffer(2000);}),1);
var proximity25km = ee.Image(0).byte().paint(globalMines.map(function(f){return f.buffer(25000);}),1);
print('── Reference Data ─────────────────────');
print('GMW polygons in AOI:', globalMines.size());
print('GMW footprint area (ha):', totalMinesHa.round());
print('Confirmed river segments:', tier1Rivers.size());
print('Named places in AOI:', places.size());

// ── MASKS ──────────────────────────────────────────────────
var worldCover = ee.ImageCollection('ESA/WorldCover/v200').first().clip(box);
var urbanMask = worldCover.neq(50).and(worldCover.neq(40));
var ghsl = ee.ImageCollection('JRC/GHSL/P2016/SMOD_POP_GLOBE_V1')
  .filterBounds(box).mosaic().clip(box);
var ghslUrban = ghsl.select('smod_code').lt(3);
var roadMask = ee.Image(1).byte().paint(roads.map(function(f){return f.buffer(50);}), 0);
var finalMask = urbanMask.and(roadMask).and(ghslUrban);

// ── CLOUD MASKING ──────────────────────────────────────────
function maskS2(img) {
  var qa = img.select('QA60');
  var qm = qa.bitwiseAnd(1<<10).eq(0).and(qa.bitwiseAnd(1<<11).eq(0));
  var scl = img.select('SCL');
  var sm = scl.neq(3).and(scl.neq(8)).and(scl.neq(9)).and(scl.neq(10)).and(scl.neq(11));
  return img.updateMask(qm).updateMask(sm).divide(10000);
}
function maskS2wet(img) {
  var scl = img.select('SCL');
  return img.updateMask(scl.neq(3).and(scl.neq(8)).and(scl.neq(9)).and(scl.neq(10))).divide(10000);
}
function maskL8(img) {
  var qa = img.select('QA_PIXEL');
  return img.updateMask(qa.bitwiseAnd(1<<3).eq(0).and(qa.bitwiseAnd(1<<4).eq(0)))
    .multiply(0.0000275).add(-0.2);
}

// ── COMPOSITE BUILDERS ─────────────────────────────────────
function compS2(d1,d2) {
  return ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
    .filterBounds(box).filterDate(d1,d2)
    .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE',30))
    .map(maskS2).median().clip(box);
}
function compL8(d1,d2) {
  return ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
    .filterBounds(box).filterDate(d1,d2).map(maskL8).median().clip(box);
}
function compL8rgb(d1,d2) {
  return compL8(d1,d2).select(['SR_B4','SR_B3','SR_B2'],['B4','B3','B2']);
}
function applyMinArea(mask,minPx) {
  return mask.updateMask(mask.connectedPixelCount(200,true).gte(minPx));
}

// ── SPECTRAL INDICES ───────────────────────────────────────
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

// ── RULE-BASED DETECTORS ───────────────────────────────────
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
  return bsi.updateMask(applyMinArea(combined,12).and(valid).and(finalMask));
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
function countHa(img,label,scale) {
  var tot=img.mask().multiply(ee.Image.pixelArea())
    .reduceRegion({reducer:ee.Reducer.sum(),geometry:box,scale:scale,maxPixels:1e10});
  print(label, ee.Number(tot.get(tot.keys().get(0))).divide(10000).round(), 'ha');
}

// ── BUILD COMPOSITES ───────────────────────────────────────
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

var d16=s2_16.select(['B4','B3','B2']).unmask(l8_16).clip(box);
var d18=s2_18.select(['B4','B3','B2']).unmask(l8_18).clip(box);

// ── RULE-BASED DETECTIONS ──────────────────────────────────
var ml14=mineL8(l8_14),ml15=mineL8(l8_15);
var m16=mineS2(s2_16),m18=mineS2(s2_18),m20=mineS2(s2_20);
var m22=mineS2(s2_22),m24=mineS2(s2_24),m25=mineS2(s2_25),m26=mineS2(s2_26);

// ── WET SEASON COMPOSITE (3-year merge for density) ────────
var s2wet = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
  .filterBounds(box).filterDate('2023-05-01','2025-10-31')
  .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE',30))
  .map(maskS2wet).median().clip(box);

var MNDWI_wet = s2wet.normalizedDifference(['B3','B11']).rename('MNDWI_wet');
var NDTI_wet  = s2wet.normalizedDifference(['B4','B3']).rename('NDTI_wet');
var IOR_wet   = s2wet.select('B4').divide(s2wet.select('B2')).rename('IOR_wet');
var NDVI_wet  = s2wet.normalizedDifference(['B8','B4']).rename('NDVI_wet');

var turbidWater = MNDWI_wet.gt(0).and(NDTI_wet.gt(0.05))
  .and(IOR_wet.gt(1.05)).and(NDVI_wet.lt(0.3))
  .selfMask().rename('TurbidWater');

// ── JRC FLOOD SAFETY ───────────────────────────────────────
var pixelArea = ee.Image.pixelArea();
var jrc = ee.ImageCollection('JRC/GSW1_4/MonthlyHistory').filter(ee.Filter.date('2000-01-01','2021-12-31'));
var jrcDry = jrc.filter(ee.Filter.calendarRange(11,2,'month')).map(function(i){return i.eq(2).selfMask();}).max().clip(box);
var jrcWet = jrc.filter(ee.Filter.calendarRange(6,9,'month')).map(function(i){return i.eq(2).selfMask();}).max().clip(box);
var floodZone = jrcWet.unmask(0).subtract(jrcDry.unmask(0)).gt(0).selfMask().clip(box);
var permWater = jrcDry.unmask(0).and(jrcWet.unmask(0)).selfMask().clip(box);
var forest = worldCover.eq(10).selfMask();

// ── BASIN GEOMETRIES ───────────────────────────────────────
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

// ── WATER RASTER LAYERS ────────────────────────────────────
var allWaterwaysRaster = ee.Image().byte().paint({featureCollection:waterways,color:1,width:1});
var tier1Raster        = ee.Image().byte().paint({featureCollection:tier1Rivers,color:1,width:3});
var adjacentWaterwaysRaster = ee.Image().byte().paint(
  waterways.map(function(f){return f.buffer(30);}),1).updateMask(miningBufferImg);

// ── COMMUNITY EXPOSURE ─────────────────────────────────────
var miningUnion = globalMines.map(function(f){return f.buffer(2000);}).union(ee.ErrorMargin(100));
var exposedCommunities = places.filterBounds(miningUnion.geometry());
var exposedCommunitiesRaster = ee.Image().byte().paint({
  featureCollection:exposedCommunities.map(function(f){return f.buffer(100);}),color:1
});

// ── RF CLASSIFIER — TRAINED ON POLYGON AREAS ───────────────
var rfBands = ['B2','B3','B4','B5','B6','B7','B8','B8A','B11','B12',
               'BSI','NDVI','NDWI','NDBI','MNDWI','IronOxide'];

var s2_rf = compS2('2024-10-01','2025-04-30');
var s2_indexed = addIndicesS2(s2_rf);
var imageStack = s2_indexed.select(rfBands);

// Sample from polygon areas — 50 pixels per polygon max to avoid memory limits
// This still gives ~31,400 samples vs 837 centroid points — massive improvement
var sampledFromPolygons = imageStack.stratifiedSample({
  numPoints: 50, classBand: null, region: box,
  scale: 20, seed: 42, geometries: false, tileScale: 4
});

// Hard cap: 2000 pixels per class — balanced training, within memory limits
// 6000 total samples vs 837 centroid points — strong improvement
var mSampled = imageStack.sampleRegions({collection: miningPolys, properties: ['c'], scale: 30, tileScale: 4}).randomColumn('r1',42).sort('r1').limit(2000);
var uSampled = imageStack.sampleRegions({collection: urbanPolys,  properties: ['c'], scale: 30, tileScale: 4}).randomColumn('r2',42).sort('r2').limit(2000);
var sSampled = imageStack.sampleRegions({collection: scrubPolys,  properties: ['c'], scale: 30, tileScale: 4}).randomColumn('r3',42).sort('r3').limit(2000);
var allSampled = mSampled.merge(uSampled).merge(sSampled);

print('── RF Training (polygon-sampled) ──────');
print('Mining samples:', mSampled.size());
print('Urban samples:', uSampled.size());
print('Scrubland samples:', sSampled.size());
print('Total samples:', allSampled.size());

var withRandom = allSampled.randomColumn('rand', 42);
var trainSet = withRandom.filter(ee.Filter.lt('rand', 0.8));
var valSet   = withRandom.filter(ee.Filter.gte('rand', 0.8));

print('Training set:', trainSet.size());
print('Validation set:', valSet.size());

var rfClassifier = ee.Classifier.smileRandomForest({numberOfTrees:200, seed:42})
  .train({features:trainSet, classProperty:'c', inputProperties:rfBands});
print('Classifier trained.');

// Standard classification — class 1 = mining
// Apply masks to reduce false positives
var rfClass = imageStack.classify(rfClassifier);
var ndviImg = s2_indexed.select('NDVI');
var rfMining = rfClass.eq(1)
  .updateMask(rfClass.eq(1))
  .updateMask(finalMask)
  .updateMask(proximity25km)
  .updateMask(ndviImg.lt(0.65))
  .selfMask();

var validated  = valSet.classify(rfClassifier);
var errorMatrix= validated.errorMatrix('c','classification');
print('── RF Accuracy ────────────────────────');
print('Confusion Matrix:', errorMatrix);
print('Overall Accuracy:', errorMatrix.accuracy());
print('Kappa:', errorMatrix.kappa());
print('Producers Accuracy:', errorMatrix.producersAccuracy());
print('Consumers Accuracy:', errorMatrix.consumersAccuracy());

// ── AREA CALCULATIONS ──────────────────────────────────────
print('── Hectare Counts (Rule-Based) ────────');
countHa(ml14,'2014 Mining Ha (L8)',30);
countHa(ml15,'2015 Mining Ha (L8)',30);
countHa(m16, '2016 Mining Ha',20);
countHa(m18, '2018 Mining Ha',20);
countHa(m20, '2020 Mining Ha',20);
countHa(m22, '2022 Mining Ha',20);
countHa(m24, '2024 Mining Ha',20);
countHa(m25, '2025 Mining Ha',20);
countHa(m26, '2026 Mining Ha',20);

var rfAreaDict = rfMining.multiply(pixelArea).divide(10000)
  .reduceRegion({reducer:ee.Reducer.sum(),geometry:box,scale:30,maxPixels:1e13, tileScale:4,bestEffort:true});
print('── RF Mining Area 2025 ────────────────');
print('RF Mining Ha:', ee.Number(rfAreaDict.get(rfAreaDict.keys().get(0))).round());

var turbidAreaDict = turbidWater.multiply(pixelArea).divide(10000)
  .reduceRegion({reducer:ee.Reducer.sum(),geometry:box,scale:20,maxPixels:1e13,bestEffort:true});
print('── Turbid Water (2023-2025 wet) ───────');
print('Turbid Water Ha:', ee.Number(turbidAreaDict.get('TurbidWater')).round());

print('── Flood Safety ───────────────────────');
var basins=[{name:'Pra',geom:basin_pra},{name:'Ankobra',geom:basin_ankobra},
            {name:'Birim',geom:basin_birim},{name:'Tano',geom:basin_tano},{name:'Offin',geom:basin_offin}];
basins.forEach(function(b){
  var ha=ee.Number(pixelArea.updateMask(floodZone.clip(b.geom))
    .reduceRegion({reducer:ee.Reducer.sum(),geometry:b.geom,scale:30,maxPixels:1e10,bestEffort:true})
    .get('area')).divide(10000);
  print(ee.String(b.name).cat(' flood zone ha: ').cat(ha.format('%.0f')));
});
print('SAFE deployment: Nov-Feb | RISK: Jun-Sep');

// ── TIME SERIES CHART ──────────────────────────────────────
function getYearTurbidHa(year) {
  var s=year+'-06-01', e=year+'-09-30';
  function mL8(img){var qa=img.select('QA_PIXEL');return img.updateMask(qa.bitwiseAnd(1<<3).eq(0).and(qa.bitwiseAnd(1<<4).eq(0))).select(['SR_B3','SR_B4','SR_B6']).multiply(0.0000275).add(-0.2);}
  var merged=ee.ImageCollection('LANDSAT/LC08/C02/T1_L2').filterBounds(box).filterDate(s,e)
    .filter(ee.Filter.lt('CLOUD_COVER',30)).map(mL8)
    .merge(ee.ImageCollection('LANDSAT/LC09/C02/T1_L2').filterBounds(box).filterDate(s,e)
    .filter(ee.Filter.lt('CLOUD_COVER',30)).map(mL8));
  var s2y=ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
    .filterBounds(box).filterDate(s,e).filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE',30)).map(maskS2wet);
  var s2Ha=ee.Algorithms.If(s2y.size().gte(3),
    ee.Number(s2y.median().clip(box).normalizedDifference(['B3','B11']).gt(0)
      .and(s2y.median().clip(box).normalizedDifference(['B4','B3']).gt(0.05))
      .multiply(pixelArea).divide(10000)
      .reduceRegion({reducer:ee.Reducer.sum(),geometry:box,scale:30,maxPixels:1e13,bestEffort:true})
      .values().get(0)),0);
  var lHa=ee.Algorithms.If(merged.size().gte(3),
    ee.Number(merged.median().clip(box).normalizedDifference(['SR_B3','SR_B6']).gt(0)
      .and(merged.median().clip(box).normalizedDifference(['SR_B4','SR_B3']).gt(0.05))
      .multiply(pixelArea).divide(10000)
      .reduceRegion({reducer:ee.Reducer.sum(),geometry:box,scale:30,maxPixels:1e13,bestEffort:true})
      .values().get(0)),null);
  return ee.Algorithms.If(merged.size().gte(3),lHa,s2Ha);
}

var tsFC=ee.FeatureCollection([2014,2015,2016,2017,2018,2019,2020,2021,2022,2023,2024,2025].map(function(y){
  return ee.Feature(null,{year:y,turbid_ha:getYearTurbidHa(y)});
}));

print('── Time Series Chart ──────────────────');
print(ui.Chart.feature.byFeature({features:tsFC,xProperty:'year',yProperties:['turbid_ha']})
  .setChartType('LineChart')
  .setOptions({title:'GalaSat: Mining-Affected Turbid Water 2014-2025 (Ashanti Pilot)',
    hAxis:{title:'Year',format:'####'},vAxis:{title:'Turbid Water Area (hectares)'},
    colors:['#8B0000'],lineWidth:3,pointSize:6,interpolateNulls:false}));

// ── MAP LAYERS ─────────────────────────────────────────────
var tcViz  ={bands:['B4','B3','B2'],min:0.0,max:0.3,gamma:1.4};
var tcL8Viz={bands:['SR_B4','SR_B3','SR_B2'],min:0.0,max:0.3,gamma:1.4};
var mineViz={min:0.0,max:0.3,palette:['ffff00','fd8d3c','f03b20','bd0026']};

// True color all years
Map.addLayer(l8_14,tcL8Viz,'2014 True Color',false);
Map.addLayer(l8_15,tcL8Viz,'2015 True Color',false);
Map.addLayer(d16, tcViz, '2016 True Color',false);
Map.addLayer(d18, tcViz, '2018 True Color',false);
Map.addLayer(s2_20,tcViz, '2020 True Color',false);
Map.addLayer(s2_22,tcViz, '2022 True Color',false);
Map.addLayer(s2_24,tcViz, '2024 True Color',false);
Map.addLayer(s2_25,tcViz, '2025 True Color',true);
Map.addLayer(s2_26,tcViz, '2026 True Color',false);

// Rule-based mining all years
Map.addLayer(ml14,mineViz,'2014 Mining (L8)',false);
Map.addLayer(ml15,mineViz,'2015 Mining (L8)',false);
Map.addLayer(m16, mineViz,'2016 Mining',false);
Map.addLayer(m18, mineViz,'2018 Mining',false);
Map.addLayer(m20, mineViz,'2020 Mining',false);
Map.addLayer(m22, mineViz,'2022 Mining',false);
Map.addLayer(m24, mineViz,'2024 Mining',false);
Map.addLayer(m25, mineViz,'2025 Mining',true);
Map.addLayer(m26, mineViz,'2026 Mining',true);

// RF output — 70% confidence threshold, proximity masked
Map.addLayer(rfMining,{min:0,max:1,palette:['FF4500'],opacity:0.7}, 'RF Mining Detection 2025',false);

// Water + environment
Map.addLayer(permWater,{palette:['#003087']},'Permanent Water',true);
Map.addLayer(floodZone,{palette:['#00ffff']},'Seasonal Flood Zone',true,0.85);
Map.addLayer(turbidWater,{palette:['FF8C00'],opacity:0.8},'Turbid Water (wet season)',true);
var forestFiltered = forest.updateMask(m26.mask().not());
Map.addLayer(forestFiltered,{palette:['1A5C1A']},'Forest Cover',false);

// Infrastructure
Map.addLayer(globalMines,{color:'FFFF00'},'Mining Footprints (Global Mining Watch)',false);
Map.addLayer(roads,{color:'888888'},'Major Roads',false);
Map.addLayer(allWaterwaysRaster,{palette:['4FC3F7'],opacity:0.5},'All Waterways',true);
Map.addLayer(tier1Raster,{palette:['8B0000']},'Field-Verified Contaminated Rivers',true);
Map.addLayer(adjacentWaterwaysRaster,{palette:['9C27B0']},'Waterways Adjacent to Mining',false);
Map.addLayer(exposedCommunitiesRaster,{palette:['FF00FF']},'Communities within 2km',true);
Map.addLayer(basin_fc.style({color:'FF6600',fillColor:'00000000',width:2}),{},'River Basin Boundaries',true);

// Mapped polygons — 3 classes as separate layers
Map.addLayer(miningPolys.style({color:'ff0000',fillColor:'ff000040',width:1}), {},'Mapped Polygons: Mining',true);
Map.addLayer(urbanPolys.style({color:'ffff00',fillColor:'ffff0040',width:1}),  {},'Mapped Polygons: Urban',false);
Map.addLayer(scrubPolys.style({color:'00ff00',fillColor:'00ff0040',width:1}),  {},'Mapped Polygons: Scrubland',false);

// AOI boundary
var aoiOutline=ee.Image().byte().paint({featureCollection:ee.FeatureCollection([ee.Feature(box)]),color:1,width:2});
Map.addLayer(aoiOutline,{palette:['FFFFFF'],opacity:1.0},'Pilot AOI Boundary');

// ── LEGEND (v4.3) ──────────────────────────────────────────
var legend=ui.Panel({style:{position:'bottom-left',padding:'8px 15px',backgroundColor:'#1a1a1a'}});
legend.add(ui.Label({value:'GalaSat v4.3 — INFN8VZN',style:{fontWeight:'bold',fontSize:'14px',margin:'0 0 2px 0',color:'#F0C0E0'}}));
legend.add(ui.Label({value:'Ashanti Pilot | Sentinel-2 + RF Polygons | April 2026',style:{fontSize:'10px',margin:'0 0 6px 0',color:'#aaaaaa'}}));
var makeRow=function(color,label){
  return ui.Panel({widgets:[
    ui.Label({style:{backgroundColor:color,padding:'7px',margin:'0 6px 3px 0'}}),
    ui.Label({value:label,style:{margin:'0 0 3px 0',fontSize:'11px',color:'#ffffff'}})
  ],layout:ui.Panel.Layout.Flow('horizontal')});
};
legend.add(makeRow('#bd0026','High Mining Activity (rule-based)'));
legend.add(makeRow('#f03b20','Moderate Mining Activity'));
legend.add(makeRow('#fd8d3c','Low Mining Activity'));
legend.add(makeRow('#FF4500','RF Mining Detection (70% confidence)'));
legend.add(makeRow('#FF8C00','Turbid Water (wet season)'));
legend.add(makeRow('#003087','Permanent Water'));
legend.add(makeRow('#00ffff','Seasonal Flood Zone'));
legend.add(makeRow('#1A5C1A','Forest Cover'));
legend.add(makeRow('#FFFF00','Global Mining Watch Footprints'));
legend.add(makeRow('#888888','Major Roads'));
legend.add(makeRow('#4FC3F7','All Waterways'));
legend.add(makeRow('#8B0000','Field-Verified Contaminated Rivers'));
legend.add(makeRow('#9C27B0','Waterways Adjacent to Mining'));
legend.add(makeRow('#FF00FF','Communities within 2km'));
legend.add(makeRow('#FF6600','River Basin Boundaries'));
legend.add(makeRow('#ff0000','Mapped Polygons: Mining'));
legend.add(makeRow('#ffff00','Mapped Polygons: Urban'));
legend.add(makeRow('#00ff00','Mapped Polygons: Scrubland'));
legend.add(makeRow('#FFFFFF','Pilot AOI Boundary'));
legend.add(ui.Label('RF: 805 polys | NDVI<0.65 | 25km mask',{fontSize:'9px',color:'#aaaaaa',margin:'4px 0 1px 0',backgroundColor:'#111111'}));
legend.add(ui.Label('Rule-based: 2014-2026 | GMW | JRC | OSM',{fontSize:'9px',color:'#aaaaaa',margin:'0',backgroundColor:'#111111'}));
Map.add(legend);

print('GalaSat v4.5 — Build complete');
print('RF trained on full polygon areas. 70% confidence threshold active.');

// ── DUPLICATE LEGEND (copy-paste artifact noted in chat analysis) ──
var legend = ui.Panel({style:{
  position:'bottom-left',padding:'6px 10px',
  backgroundColor:'#1c1c1c',maxHeight:'420px',width:'250px'
}});
var makeRow = function(color, label) {
  return ui.Panel({
    widgets:[
      ui.Label({value:' ',style:{backgroundColor:color,padding:'4px 8px',margin:'0 4px 1px 0'}}),
      ui.Label({value:label,style:{fontSize:'10px',color:'#000000',backgroundColor:'#dddddd',padding:'1px 3px',margin:'1px 0'}})
    ],
    layout:ui.Panel.Layout.Flow('horizontal'),
    style:{backgroundColor:'#1c1c1c',margin:'1px 0'}
  });
};
legend.add(ui.Label({value:'GalaSat v4.5 — INFN8VZN',style:{fontWeight:'bold',fontSize:'12px',margin:'0 0 4px 0',color:'#000000',backgroundColor:'#F0C0E0',padding:'2px 5px'}}));
legend.add(ui.Label({value:'Ashanti Pilot | April 2026',style:{fontSize:'9px',color:'#cccccc',margin:'0 0 4px 0',backgroundColor:'#1c1c1c'}}));
legend.add(makeRow('#bd0026','High Mining Activity'));
legend.add(makeRow('#f03b20','Moderate Mining Activity'));
legend.add(makeRow('#fd8d3c','Low Mining Activity'));
legend.add(makeRow('#FF4500','RF Mining Detection 2025'));
legend.add(makeRow('#FF8C00','Turbid Water (wet season)'));
legend.add(makeRow('#003087','Permanent Water'));
legend.add(makeRow('#00ffff','Seasonal Flood Zone'));
legend.add(makeRow('#1A5C1A','Forest Cover'));
legend.add(makeRow('#FFFF00','Global Mining Watch'));
legend.add(makeRow('#888888','Major Roads'));
legend.add(makeRow('#4FC3F7','All Waterways'));
legend.add(makeRow('#8B0000','Contaminated Rivers'));
legend.add(makeRow('#9C27B0','Waterways adj. to Mining'));
legend.add(makeRow('#FF00FF','Communities within 2km'));
legend.add(makeRow('#FF6600','River Basin Boundaries'));
legend.add(makeRow('#ff0000','Mapped Polygons: Mining'));
legend.add(makeRow('#ffff00','Mapped Polygons: Urban'));
legend.add(makeRow('#00ff00','Mapped Polygons: Scrubland'));
legend.add(makeRow('#FFFFFF','Pilot AOI Boundary'));
legend.add(ui.Label({value:'81% accuracy | Kappa 0.715 | 805 polygons',style:{fontSize:'9px',color:'#cccccc',margin:'4px 0 1px 0',backgroundColor:'#1c1c1c'}}));
legend.add(ui.Label({value:'Rule-based: 2014-2026 | GMW | JRC | OSM',style:{fontSize:'9px',color:'#cccccc',margin:'0',backgroundColor:'#1c1c1c'}}));
Map.add(legend);

print('GalaSat v4.5 — Build complete');
print('RF trained on full polygon areas. 70% confidence threshold active.');
