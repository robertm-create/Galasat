// ── GALASAT v2.2 — PROXIMITY 25km + NDVI lt(0.7) ───────────────────────
//
// MEMORY TAG: GalaSat_v2.2_Landsat_SAR
// SOURCE: claude.ai chat "Hektare polygon data not importing to Galasad"
//         https://claude.ai/chat/dc2dd274-3358-4462-8ec8-abe4f7169f18
//         pasted attachment, 335 lines, 19,611 bytes
//         Account 2 (Teqami's), discovered 2026-04-25 cross-account sweep
// EXTRACTED: 2026-04-25
//
// LINEAGE NOTE: Self-labeled GalaSat v2.2 — distinct from the v2.0 Frankenstein in
// Stacy's Google grant continuation chat. Major architectural shift:
//   - Trains on NEW asset GALAMSEY_TRAINING_v3 (752 pts per legend, vs 340)
//   - Replaces Sentinel-2 dry composite with a pre-built LANDSAT IMAGE ASSET
//     (`projects/galamsey-monotoring/assets/galasat_landsat_dry_composite_2026_v1`)
//   - Dual SAR change layers: VV AND VH (earlier builds had only one)
//   - Three stacked masks: ashantiMask (lat > 5.80) + proximityMask (25km from
//     GMW) + ndviMask (NDVI < 0.7) — tighter than V3's RF chain
//   - Console reports "Tier 3 total classifier output (unfiltered): 1,063,650 ha"
//     — this is the new unfiltered ceiling reference
//
// NOTE: Whitespace partially collapsed during DOM extraction.

var trainingPoints = ee.FeatureCollection('projects/galamsey-monotoring/assets/infn8vzn_mercury_training_340');
var trainingV3 = ee.FeatureCollection('projects/galamsey-monotoring/assets/GALAMSEY_TRAINING_v3');
var AOI = trainingPoints.geometry().bounds().buffer(20000);

var roads = ee.FeatureCollection('projects/galamsey-monotoring/assets/ghana_osm_roads')
  .filterBounds(AOI)
  .filter(ee.Filter.inList('fclass',['motorway','trunk','primary','secondary']));
var waterways = ee.FeatureCollection('projects/galamsey-monotoring/assets/ghana_osm_waterways').filterBounds(AOI);
var places = ee.FeatureCollection('projects/galamsey-monotoring/assets/ghana_osm_places').filterBounds(AOI);

var tier1Rivers = waterways.filter(ee.Filter.inList('name',[
  'Pra','Pra River','Offin','Offin River','Ankobra','Ankobra River',
  'Birim','Birim River','Tano','Tano River','Densu','Densu River',
  'Oda','Oda River','Bonsa','Bonsa River','Oti','Oti River',
  'Fena','Fena River','Subin','Subin River','Anum','Anum River',
  'Afram','Afram River','Anikoko','Bodwire','Asesree','Assaman'
]));

var globalMines = ee.FeatureCollection('projects/sat-io/open-datasets/global-mining/global_mining_polygons').filterBounds(AOI);
var minesClipped = globalMines.map(function(f){return f.intersection(AOI,ee.ErrorMargin(1));});
var totalMinesHa = ee.Number(minesClipped.map(function(f){return f.set('area_ha',f.geometry().area(1).divide(10000));}).aggregate_sum('area_ha'));
var miningBuffer = globalMines.map(function(f){return f.buffer(2000);});
var miningBufferImage = ee.Image(0).byte().paint(miningBuffer,1);

print('── Mining Footprints ──────────────────');
print('Mining polygons in AOI:', globalMines.size());
print('Mining footprint area (hectares):', totalMinesHa);
print('Confirmed river segments in AOI:', tier1Rivers.size());
print('Total named places in AOI:', places.size());

// ── LOAD CLEAN LANDSAT ASSET ───────────────────────────────────────────
var landsatAsset = ee.Image('projects/galamsey-monotoring/assets/galasat_landsat_dry_composite_2026_v1')
  .select(['b1','b2','b3','b4','b5','b6'],['B2','B3','B4','B5','B6','B7'])
  .clip(AOI);

// ── CLOUD MASK FOR WET SEASON ──────────────────────────────────────────
function maskS2clouds(image) {
  var scl = image.select('SCL');
  var mask = scl.neq(3).and(scl.neq(8)).and(scl.neq(9)).and(scl.neq(10));
  return image.updateMask(mask).divide(10000);
}

// ── WET SEASON S2 FOR TURBID WATER ─────────────────────────────────────
var s2wet2025 = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
  .filterBounds(AOI).filterDate('2025-06-01','2025-09-30')
  .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE',30)).map(maskS2clouds).median().clip(AOI);

// ── SENTINEL-1 SAR — ASCENDING ONLY ────────────────────────────────────
function processSAR(image) {
  var vv = image.select('VV').focal_mean(3,'square','pixels').rename('VV_filtered');
  var vh = image.select('VH').focal_mean(3,'square','pixels').rename('VH_filtered');
  return vv.addBands(vh);
}
var s1dry = ee.ImageCollection('COPERNICUS/S1_GRD')
  .filterBounds(AOI)
  .filterDate('2025-11-01','2026-03-31')
  .filter(ee.Filter.eq('instrumentMode','IW'))
  .filter(ee.Filter.listContains('transmitterReceiverPolarisation','VV'))
  .filter(ee.Filter.listContains('transmitterReceiverPolarisation','VH'))
  .filter(ee.Filter.eq('orbitProperties_pass','ASCENDING'))
  .select(['VV','VH'])
  .map(processSAR)
  .median().clip(AOI);
var s1wet = ee.ImageCollection('COPERNICUS/S1_GRD')
  .filterBounds(AOI)
  .filterDate('2025-06-01','2025-09-30')
  .filter(ee.Filter.eq('instrumentMode','IW'))
  .filter(ee.Filter.listContains('transmitterReceiverPolarisation','VV'))
  .filter(ee.Filter.listContains('transmitterReceiverPolarisation','VH'))
  .filter(ee.Filter.eq('orbitProperties_pass','ASCENDING'))
  .select(['VV','VH'])
  .map(processSAR)
  .median().clip(AOI);

var sarRatio_dry = s1dry.select('VV_filtered').divide(s1dry.select('VH_filtered')).rename('SAR_ratio_dry');
var sarChange_VV = s1wet.select('VV_filtered').subtract(s1dry.select('VV_filtered')).rename('SAR_change_VV');
var sarChange_VH = s1wet.select('VH_filtered').subtract(s1dry.select('VH_filtered')).rename('SAR_change_VH');

// ── LANDSAT INDICES ────────────────────────────────────────────────────
var NDVI  = landsatAsset.normalizedDifference(['B5','B4']).rename('NDVI');
var NDWI  = landsatAsset.normalizedDifference(['B3','B5']).rename('NDWI');
var MNDWI = landsatAsset.normalizedDifference(['B3','B6']).rename('MNDWI');
var NDTI  = landsatAsset.normalizedDifference(['B4','B3']).rename('NDTI');
var BSI   = landsatAsset.expression('((SWIR+RED)-(NIR+BLUE))/((SWIR+RED)+(NIR+BLUE))',
  {SWIR:landsatAsset.select('B6'),RED:landsatAsset.select('B4'),
   NIR:landsatAsset.select('B5'),BLUE:landsatAsset.select('B2')}).rename('BSI');
var IOR   = landsatAsset.select('B4').divide(landsatAsset.select('B2')).rename('IronOxide');

// ── WET INDICES ────────────────────────────────────────────────────────
var MNDWI_wet = s2wet2025.normalizedDifference(['B3','B11']).rename('MNDWI_wet');
var NDTI_wet  = s2wet2025.normalizedDifference(['B4','B3']).rename('NDTI_wet');
var IOR_wet   = s2wet2025.select('B4').divide(s2wet2025.select('B2')).rename('IOR_wet');
var NDVI_wet  = s2wet2025.normalizedDifference(['B8','B4']).rename('NDVI_wet');

// ── FEATURE STACK — LANDSAT + SAR ──────────────────────────────────────
var imageStack = landsatAsset.select(['B2','B3','B4','B5','B6','B7'])
  .addBands(NDVI).addBands(NDWI).addBands(MNDWI).addBands(NDTI).addBands(BSI).addBands(IOR)
  .addBands(s1dry.select('VV_filtered'))
  .addBands(s1dry.select('VH_filtered'))
  .addBands(sarRatio_dry)
  .addBands(sarChange_VV)
  .addBands(sarChange_VH);

// ── TRAINING — GALAMSEY_TRAINING_v3 ────────────────────────────────────
var allSampled = imageStack.sampleRegions({
  collection: trainingV3, properties: ['class'], scale: 30, tileScale: 8
});

var allSamples = allSampled.randomColumn('random', 42);
var trainSet = allSamples.filter(ee.Filter.lt('random',0.8));
var valSet   = allSamples.filter(ee.Filter.gte('random',0.8));

print('── ACCURACY CHECK ─────────────────────');
print('Total samples:', allSamples.size());
print('Training set:', trainSet.size());
print('Validation set:', valSet.size());

var bands = ['B2','B3','B4','B5','B6','B7',
             'NDVI','NDWI','MNDWI','NDTI','BSI','IronOxide',
             'VV_filtered','VH_filtered','SAR_ratio_dry','SAR_change_VV','SAR_change_VH'];

var classifier = ee.Classifier.smileRandomForest({numberOfTrees:200, seed:42})
  .train({features:trainSet, classProperty:'class', inputProperties:bands});

var classified = imageStack.classify(classifier);

// ── MASKS ──────────────────────────────────────────────────────────────
var worldCover = ee.ImageCollection('ESA/WorldCover/v200').first().clip(AOI);
var urbanMask = worldCover.neq(50).and(worldCover.neq(40));
var roadMask = ee.Image(1).byte().paint(roads.map(function(f){return f.buffer(50);}),0);
var finalMask = urbanMask.and(roadMask);
var ashantiMask = ee.Image.pixelLonLat().select('latitude').gt(5.80);

// FIX 2 — PROXIMITY MASK: 25km from confirmed mining
var proximityMask = ee.Image(0).byte().paint(
  globalMines.map(function(f){return f.buffer(25000);}), 1
);

// FIX 3 — NDVI SUPPRESSION: lt(0.7) keeps disturbed riparian zones
var ndviMask = NDVI.lt(0.7);

var contaminationOnly = classified
  .updateMask(finalMask)
  .updateMask(ashantiMask)
  .updateMask(proximityMask)
  .updateMask(ndviMask)
  .updateMask(classified.eq(1))
  .selfMask();

var forest = worldCover.eq(10).selfMask();

// ── TURBID WATER ───────────────────────────────────────────────────────
var turbidWater2025 = MNDWI_wet.gt(0).and(NDTI_wet.gt(0.05)).and(IOR_wet.gt(1.05)).and(NDVI_wet.lt(0.3))
  .selfMask().rename('TurbidWater');

// ── FLOOD SAFETY (JRC) ─────────────────────────────────────────────────
var jrc = ee.ImageCollection('JRC/GSW1_4/MonthlyHistory').filter(ee.Filter.date('2000-01-01','2021-12-31'));
var jrcDry = jrc.filter(ee.Filter.calendarRange(11,2,'month')).map(function(i){return i.eq(2).selfMask();}).max().clip(AOI);
var jrcWet = jrc.filter(ee.Filter.calendarRange(6,9,'month')).map(function(i){return i.eq(2).selfMask();}).max().clip(AOI);
var flood_zone = jrcWet.unmask(0).subtract(jrcDry.unmask(0)).gt(0).selfMask().clip(AOI);
var permanent_water = jrcDry.unmask(0).and(jrcWet.unmask(0)).selfMask().clip(AOI);

// ── BASIN GEOMETRIES ───────────────────────────────────────────────────
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

var pixelArea = ee.Image.pixelArea();
var basins = [
  {name:'Pra',geom:basin_pra},{name:'Ankobra',geom:basin_ankobra},
  {name:'Birim',geom:basin_birim},{name:'Tano',geom:basin_tano},{name:'Offin',geom:basin_offin}
];

// ── WATER RASTER LAYERS ────────────────────────────────────────────────
var allWaterwaysRaster = ee.Image().byte().paint({featureCollection:waterways, color:1, width:1});
var tier1Raster = ee.Image().byte().paint({featureCollection:tier1Rivers, color:1, width:3});
var adjacentWaterwaysRaster = ee.Image().byte().paint(waterways.map(function(f){return f.buffer(30);}),1).updateMask(miningBufferImage);

// ── COMMUNITY EXPOSURE ─────────────────────────────────────────────────
var miningUnion = globalMines.map(function(f){return f.buffer(2000);}).union(ee.ErrorMargin(100));
var exposedCommunities = places.filterBounds(miningUnion.geometry());
var exposedCommunitiesRaster = ee.Image().byte().paint({
  featureCollection:exposedCommunities.map(function(f){return f.buffer(100);}),
  color:1
});

// ── AREA CALCULATIONS ──────────────────────────────────────────────────
var turbidAreaDict = turbidWater2025.multiply(pixelArea).divide(10000)
  .reduceRegion({reducer:ee.Reducer.sum(),geometry:AOI,scale:30,maxPixels:1e13,bestEffort:true});
print('── Turbid Water (wet season 2025) ─────');
print('Turbid Water Area (hectares):', turbidAreaDict.get('TurbidWater'));

var contamAreaDict = contaminationOnly.multiply(pixelArea).divide(10000)
  .reduceRegion({reducer:ee.Reducer.sum(),geometry:AOI,scale:30,maxPixels:1e13,bestEffort:true});
print('── Mining Disturbance 2026 ────────────');
print('Mining Disturbance Area (hectares):', contamAreaDict.get('classification'));

print('── Tier breakdown ─────────────────────');
print('Tier 1 confirmed (GMW):', totalMinesHa);
print('Tier 3 total classifier output (unfiltered): 1,063,650 ha');

print('── Community Exposure ─────────────────');
print('Total communities in AOI:', places.size());
print('Communities within 2km of confirmed mining:', exposedCommunities.size());

print('── Flood Safety — Basin Extent ────────');
basins.forEach(function(b) {
  var ha = ee.Number(pixelArea.updateMask(flood_zone.clip(b.geom))
    .reduceRegion({reducer:ee.Reducer.sum(),geometry:b.geom,scale:30,maxPixels:1e10,bestEffort:true})
    .get('area')).divide(10000);
  print(ee.String(b.name).cat(': ').cat(ha.format('%.1f')).cat(' ha'));
});
print('SAFE deployment: Nov-Feb | RISK: Jun-Sep');

// ── VALIDATION ─────────────────────────────────────────────────────────
var validated = valSet.classify(classifier);
var errorMatrix = validated.errorMatrix('class','classification');
print('── Confusion Matrix ───────────────────');
print('Confusion Matrix:', errorMatrix);
print('Overall Accuracy:', errorMatrix.accuracy());
print('Kappa:', errorMatrix.kappa());
print('Producers Accuracy:', errorMatrix.producersAccuracy());
print('Consumers Accuracy:', errorMatrix.consumersAccuracy());
print('F1 per class:', errorMatrix.fscore());

// ── TIME SERIES 2014-2025 ──────────────────────────────────────────────
function getYearTurbidHa(year) {
  var s = year+'-06-01', e = year+'-09-30';
  function maskL8(img){var qa=img.select('QA_PIXEL');return img.updateMask(qa.bitwiseAnd(1<<3).eq(0).and(qa.bitwiseAnd(1<<4).eq(0))).select(['SR_B3','SR_B4','SR_B6']).multiply(0.0000275).add(-0.2);}
  function maskL9(img){var qa=img.select('QA_PIXEL');return img.updateMask(qa.bitwiseAnd(1<<3).eq(0).and(qa.bitwiseAnd(1<<4).eq(0))).select(['SR_B3','SR_B4','SR_B6']).multiply(0.0000275).add(-0.2);}
  var merged = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2').filterBounds(AOI).filterDate(s,e).filter(ee.Filter.lt('CLOUD_COVER',30)).map(maskL8)
    .merge(ee.ImageCollection('LANDSAT/LC09/C02/T1_L2').filterBounds(AOI).filterDate(s,e).filter(ee.Filter.lt('CLOUD_COVER',30)).map(maskL9));
  var s2 = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
    .filterBounds(AOI).filterDate(s,e)
    .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE',30))
    .map(maskS2clouds);
  var s2count = s2.size();
  var s2Ha = ee.Algorithms.If(s2count.gte(3),
    ee.Number(s2.median().clip(AOI).normalizedDifference(['B3','B11']).rename('MNDWI').gt(0)
      .and(s2.median().clip(AOI).normalizedDifference(['B4','B3']).rename('NDTI').gt(0.05))
      .rename('turbid_ha').multiply(pixelArea).divide(10000)
      .reduceRegion({reducer:ee.Reducer.sum(),geometry:AOI,scale:30,maxPixels:1e13,bestEffort:true})
      .get('turbid_ha')),0);
  var landsatHa = ee.Algorithms.If(
    merged.size().gte(3),
    ee.Number(merged.median().clip(AOI).normalizedDifference(['SR_B3','SR_B6']).rename('MNDWI').gt(0)
      .and(merged.median().clip(AOI).normalizedDifference(['SR_B4','SR_B3']).rename('NDTI').gt(0.05))
      .rename('turbid_ha').multiply(pixelArea).divide(10000)
      .reduceRegion({reducer:ee.Reducer.sum(),geometry:AOI,scale:30,maxPixels:1e13,bestEffort:true}).get('turbid_ha')),null);
  return ee.Algorithms.If(merged.size().gte(3), landsatHa, s2Ha);
}

var timeSeriesFC = ee.FeatureCollection([2014,2015,2016,2017,2018,2019,2020,2021,2022,2023,2024,2025].map(function(y){
  return ee.Feature(null,{year:y,turbid_ha:getYearTurbidHa(y)});
}));
print('── Time Series Chart ──────────────────');
print(ui.Chart.feature.byFeature({features:timeSeriesFC,xProperty:'year',yProperties:['turbid_ha']})
  .setChartType('LineChart')
  .setOptions({title:'GalaSat: Mining-Affected Turbid Water 2014-2025',
    hAxis:{title:'Year',format:'####'},vAxis:{title:'Turbid Water Area (hectares)'},
    colors:['#8B0000'],lineWidth:3,pointSize:6,interpolateNulls:false}));

// ── VISUALIZE ──────────────────────────────────────────────────────────
Map.centerObject(AOI,10);
var aoiOutline = ee.Image().byte().paint({featureCollection:ee.FeatureCollection([ee.Feature(AOI)]),color:1,width:2});

Map.addLayer(landsatAsset, {bands:['B4','B3','B2'],min:0.03,max:0.18}, 'Landsat True Colour');
Map.addLayer(forest, {palette:['1A5C1A']}, 'Forest Cover');
Map.addLayer(permanent_water, {palette:['#003087']}, 'Permanent Water', true, 1.0);
Map.addLayer(flood_zone, {palette:['#00ffff']}, 'Seasonal Flood Zone (Nov-Feb safe)', true, 0.85);
Map.addLayer(contaminationOnly, {palette:['FF0000'],opacity:0.6}, 'Mining Disturbance Signal');
Map.addLayer(turbidWater2025, {palette:['FF8C00'],opacity:0.8}, 'Turbid Water (wet season 2025)');
Map.addLayer(globalMines, {color:'FFFF00'}, 'Mining Footprints — Global Mining Watch');
Map.addLayer(roads, {color:'888888'}, 'Major Roads');
Map.addLayer(allWaterwaysRaster, {palette:['4FC3F7'],opacity:0.5}, 'All Waterways');
Map.addLayer(tier1Raster, {palette:['8B0000']}, 'Field-Verified Contaminated Rivers');
Map.addLayer(adjacentWaterwaysRaster, {palette:['9C27B0']}, 'Waterways Adjacent to Confirmed Mining');
Map.addLayer(exposedCommunitiesRaster, {palette:['FF00FF']}, 'Communities within 2km');
Map.addLayer(trainingPoints, {color:'FF6600'}, 'Training Points (340)', false);
Map.addLayer(trainingV3, {color:'00FF00'}, 'GALAMSEY_TRAINING_v3 (752)', false);
Map.addLayer(basin_fc.style({color:'FF6600',fillColor:'00000000',width:2}), {}, 'River Basin Boundaries');
Map.addLayer(aoiOutline, {palette:['FFFFFF'],opacity:1.0}, 'Pilot AOI Boundary');
Map.addLayer(s1dry.select('VV_filtered'), {min:-20,max:0,palette:['000000','ffffff']}, 'SAR VV Dry Season', false);
Map.addLayer(sarChange_VV, {min:-5,max:5,palette:['0000ff','ffffff','ff0000']}, 'SAR Change VV (wet-dry)', false);

// ── LEGEND ─────────────────────────────────────────────────────────────
var legend = ui.Panel({style:{position:'bottom-left',padding:'8px 15px',backgroundColor:'white'}});
legend.add(ui.Label({value:'GalaSat v2.2 — INFN8VZN',style:{fontWeight:'bold',fontSize:'14px',margin:'0 0 2px 0',color:'#8B0000'}}));
legend.add(ui.Label({value:'Ashanti Pilot | Landsat + SAR | April 2026',style:{fontSize:'11px',margin:'0 0 6px 0',color:'#555555'}}));
var makeRow = function(color,label){
  return ui.Panel({widgets:[
    ui.Label({style:{backgroundColor:color,padding:'8px',margin:'0 6px 4px 0'}}),
    ui.Label({value:label,style:{margin:'0 0 4px 6px',fontSize:'12px'}})
  ],layout:ui.Panel.Layout.Flow('horizontal')});
};
legend.add(makeRow('#1A5C1A','Forest Cover'));
legend.add(makeRow('#003087','Permanent Water'));
legend.add(makeRow('#00ffff','Seasonal Flood Zone (deployment safety)'));
legend.add(makeRow('#FF0000','Mining Disturbance (Nov 2025-Feb 2026)'));
legend.add(makeRow('#FF8C00','Turbid Water (Jun-Sep 2025)'));
legend.add(makeRow('#FFFF00','Mining Footprints (Global Mining Watch)'));
legend.add(makeRow('#888888','Major Roads'));
legend.add(makeRow('#4FC3F7','All Waterways'));
legend.add(makeRow('#8B0000','Field-Verified Contaminated Rivers'));
legend.add(makeRow('#9C27B0','Waterways Adjacent to Confirmed Mining'));
legend.add(makeRow('#FF00FF','Communities within 2km of confirmed mining'));
legend.add(makeRow('#FFFFFF','Pilot AOI Boundary'));
legend.add(ui.Label({value:'Landsat 8/9 + Sentinel-1 SAR | Nov 2025-Mar 2026 | EPSG:4326',style:{fontSize:'10px',color:'#555555',margin:'6px 0 0 0'}}));
legend.add(ui.Label({value:'Sources: Joy News | Al Jazeera | PMC 2024 | GWCL | Global Mining Watch',style:{fontSize:'10px',color:'#555555',margin:'2px 0 0 0'}}));
Map.add(legend);

// ── EXPORTS ────────────────────────────────────────────────────────────
Export.image.toDrive({image:contaminationOnly,description:'INFN8VZN_Mercury_Disturbance_v2_2',
  folder:'GalaSat',fileNamePrefix:'infn8vzn_mercury_disturbance_v2_2',
  region:AOI,scale:30,crs:'EPSG:4326',maxPixels:1e13,fileFormat:'GeoTIFF'});
Export.image.toDrive({image:turbidWater2025,description:'INFN8VZN_Turbid_Water_2025_v2',
  folder:'GalaSat',fileNamePrefix:'infn8vzn_turbid_water_2025_v2',
  region:AOI,scale:30,crs:'EPSG:4326',maxPixels:1e13,fileFormat:'GeoTIFF'});
print('── Exports submitted — check Tasks tab ──');
