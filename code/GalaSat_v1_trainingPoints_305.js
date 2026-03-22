// ============================================================
// GalaSat v1.0 — trainingPoints architecture, 305-line build
// (Account 2 / Teqami's Hektare polygon import chat)
//
// MEMORY TAG: GalaSat_v1.0_trainingPoints_305
// SOURCE: claude.ai chat "Hektare polygon data not importing to Galasad"
//         https://claude.ai/chat/dc2dd274-3358-4462-8ec8-abe4f7169f18
//         pasted attachment, 305 lines, 18,240 bytes
// EXTRACTED: 2026-04-25 (cross-account sweep)
//
// LINEAGE NOTE: trainingPoints branch member, sized between the 286 and 311 variants
// from Stacy's chats. Adds vs the 274/259 baseline:
//   - WorldPop 2020 population density layer + popInContamZone reduction
//   - ESA WorldCover crop zones (class 40)
//   - IOR > 1.5 mercury-proxy visualization layer
//   - Larger Tano basin geometry [-3.20, 6.50, -2.00, 7.80] (matches 506 variant)
// Otherwise functionally equivalent to the trainingPoints lineage.
//
// NOTE: Whitespace partially collapsed during DOM extraction.

var trainingPoints = ee.FeatureCollection('projects/galamsey-monotoring/assets/infn8vzn_mercury_training_340');
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

// ── CLOUD MASK ─────────────────────────────────────────────────────────
function maskS2clouds(image) {
  var scl = image.select('SCL');
  var mask = scl.neq(3).and(scl.neq(8)).and(scl.neq(9)).and(scl.neq(10));
  return image.updateMask(mask).divide(10000);
}

// ── COMPOSITES ─────────────────────────────────────────────────────────
var s2dry = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
  .filterBounds(AOI).filterDate('2026-01-01','2026-02-28')
  .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE',10)).map(maskS2clouds).median().clip(AOI);
var s2wet2025 = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
  .filterBounds(AOI).filterDate('2025-06-01','2025-09-30')
  .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE',30)).map(maskS2clouds).median().clip(AOI);

// ── DRY INDICES ────────────────────────────────────────────────────────
var NDVI = s2dry.normalizedDifference(['B8','B4']).rename('NDVI');
var NDWI = s2dry.normalizedDifference(['B3','B8']).rename('NDWI');
var MNDWI = s2dry.normalizedDifference(['B3','B11']).rename('MNDWI');
var NDTI = s2dry.normalizedDifference(['B4','B3']).rename('NDTI');
var BSI = s2dry.expression('((SWIR+RED)-(NIR+BLUE))/((SWIR+RED)+(NIR+BLUE))',
  {SWIR:s2dry.select('B11'),RED:s2dry.select('B4'),NIR:s2dry.select('B8'),BLUE:s2dry.select('B2')}).rename('BSI');
var IOR = s2dry.select('B4').divide(s2dry.select('B2')).rename('IronOxide');

// ── WET INDICES ────────────────────────────────────────────────────────
var MNDWI_wet = s2wet2025.normalizedDifference(['B3','B11']).rename('MNDWI_wet');
var NDTI_wet = s2wet2025.normalizedDifference(['B4','B3']).rename('NDTI_wet');
var IOR_wet = s2wet2025.select('B4').divide(s2wet2025.select('B2')).rename('IOR_wet');
var NDVI_wet = s2wet2025.normalizedDifference(['B8','B4']).rename('NDVI_wet');

// ── FEATURE STACK ──────────────────────────────────────────────────────
var imageStack = s2dry.select(['B2','B3','B4','B5','B6','B7','B8','B8A','B11','B12'])
  .addBands(NDVI).addBands(NDWI).addBands(MNDWI).addBands(NDTI).addBands(BSI).addBands(IOR);

// ── TRAINING ───────────────────────────────────────────────────────────
var cleanPoints = imageStack.updateMask(NDVI.gt(0.6)).sample({
  region:AOI, scale:30, numPixels:278, seed:42, geometries:true
}).map(function(f){return f.set('class',0);});

var mercurySampled = imageStack.sampleRegions({
  collection:trainingPoints.filter(ee.Filter.eq('class',1)),
  properties:['class'], scale:10, tileScale:8
});

var allSamples = mercurySampled.merge(cleanPoints).randomColumn('random',42);
var trainSet = allSamples.filter(ee.Filter.lt('random',0.8));
var valSet = allSamples.filter(ee.Filter.gte('random',0.8));

print('── ACCURACY CHECK ─────────────────────');
print('Total samples:', allSamples.size());
print('Training set:', trainSet.size());
print('Validation set:', valSet.size());

var bands = ['B2','B3','B4','B5','B6','B7','B8','B8A','B11','B12','NDVI','NDWI','MNDWI','NDTI','BSI','IronOxide'];
var classifier = ee.Classifier.smileRandomForest({numberOfTrees:200, seed:42})
  .train({features:trainSet, classProperty:'class', inputProperties:bands});

var classified = imageStack.classify(classifier);

// ── MASKS ──────────────────────────────────────────────────────────────
var worldCover = ee.ImageCollection('ESA/WorldCover/v200').first().clip(AOI);
var urbanMask = worldCover.neq(50).and(worldCover.neq(40));
var roadMask = ee.Image(1).byte().paint(roads.map(function(f){return f.buffer(50);}),0);
var finalMask = urbanMask.and(roadMask);

var contaminationOnly = classified.updateMask(finalMask).updateMask(classified.eq(1)).selfMask();
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

// ── BASIN GEOMETRIES (Tano basin = larger geometry) ────────────────────
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

// ── COMMUNITY EXPOSURE — GMW polygon proximity ─────────────────────────
var miningUnion = globalMines.map(function(f){return f.buffer(2000);}).union(ee.ErrorMargin(100));
var exposedCommunities = places.filterBounds(miningUnion.geometry());
var exposedCommunitiesRaster = ee.Image().byte().paint({
  featureCollection:exposedCommunities.map(function(f){return f.buffer(100);}),
  color:1
});

// ── POPULATION DENSITY — WorldPop 2020 ─────────────────────────────────
var worldPop = ee.ImageCollection('WorldPop/GP/100m/pop')
  .filter(ee.Filter.eq('country','GHA'))
  .filter(ee.Filter.eq('year',2020))
  .first().clip(AOI);

// ── CROP ZONES — ESA WorldCover ────────────────────────────────────────
var cropZones = worldCover.eq(40).selfMask();

// ── IOR VISUALIZATION — already in feature stack ───────────────────────
var iorViz = IOR.updateMask(IOR.gt(1.5)).selfMask();

var turbidAreaDict = turbidWater2025.multiply(pixelArea).divide(10000)
  .reduceRegion({reducer:ee.Reducer.sum(),geometry:AOI,scale:30,maxPixels:1e13});
print('── Turbid Water (wet season 2025) ─────');
print('Turbid Water Area (hectares):', turbidAreaDict.get('TurbidWater'));

var contamAreaDict = contaminationOnly.multiply(pixelArea).divide(10000)
  .reduceRegion({reducer:ee.Reducer.sum(),geometry:AOI,scale:30,maxPixels:1e13});
print('── Mining Disturbance 2026 ────────────');
print('Mining Disturbance Area (hectares):', contamAreaDict.get('classification'));

// Population exposure
var popInContamZone = worldPop.updateMask(contaminationOnly);
var popDict = popInContamZone.reduceRegion({
  reducer:ee.Reducer.sum(), geometry:AOI, scale:100, maxPixels:1e13, bestEffort:true
});
print('── Population Exposure ────────────────');
print('Est. population in contamination zone:', popDict.get('population'));
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
  // Sentinel-2 fallback for years with insufficient Landsat coverage
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
      .get('turbid_ha')),
    0
  );
  var landsatHa = ee.Algorithms.If(
    merged.size().gte(3),
    ee.Number(merged.median().clip(AOI).normalizedDifference(['SR_B3','SR_B6']).rename('MNDWI').gt(0)
      .and(merged.median().clip(AOI).normalizedDifference(['SR_B4','SR_B3']).rename('NDTI').gt(0.05))
      .rename('turbid_ha').multiply(pixelArea).divide(10000)
      .reduceRegion({reducer:ee.Reducer.sum(),geometry:AOI,scale:30,maxPixels:1e13,bestEffort:true}).get('turbid_ha')),
    null
  );
  // Use Landsat if available, fall back to Sentinel-2
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
Map.setOptions('SATELLITE');
var aoiOutline = ee.Image().byte().paint({featureCollection:ee.FeatureCollection([ee.Feature(AOI)]),color:1,width:2});

// NOTE: original source references s2dry but here it should map to a true-colour layer.
// Restored as in source (preserved as captured):
Map.addLayer(s2dry, {bands:['B4','B3','B2'],min:0,max:0.3}, 'True Colour 2025-26');
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
Map.addLayer(cropZones, {palette:['#FFA500']}, 'Crop Zones (farmland)', false);
Map.addLayer(iorViz, {palette:['#8B4513'], min:1.5, max:3.0}, 'Iron Oxide Ratio (mercury proxy)', false);
Map.addLayer(worldPop, {min:0, max:200, palette:['#FFFFE0','#FFA500','#FF0000']}, 'Population Density', false);
Map.addLayer(basin_fc.style({color:'FF6600',fillColor:'00000000',width:2}), {}, 'River Basin Boundaries');
Map.addLayer(aoiOutline, {palette:['FFFFFF'],opacity:1.0}, 'Pilot AOI Boundary');

// ── LEGEND ─────────────────────────────────────────────────────────────
var legend = ui.Panel({style:{position:'bottom-left',padding:'8px 15px',backgroundColor:'white'}});
legend.add(ui.Label({value:'GalaSat v1.0 — INFN8 VZN',style:{fontWeight:'bold',fontSize:'14px',margin:'0 0 2px 0',color:'#8B0000'}}));
legend.add(ui.Label({value:'Ashanti Pilot, Ghana | March 2026',style:{fontSize:'11px',margin:'0 0 6px 0',color:'#555555'}}));
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
legend.add(makeRow('#FF6600','Training Points (340)'));
legend.add(makeRow('#FFFFFF','Pilot AOI Boundary'));
legend.add(ui.Label({value:'Sentinel-2 SR | Dual composite 2025-26 | EPSG:4326',style:{fontSize:'10px',color:'#555555',margin:'6px 0 0 0'}}));
legend.add(ui.Label({value:'Sources: Joy News | Al Jazeera | PMC 2024 | GWCL | Global Mining Watch',style:{fontSize:'10px',color:'#555555',margin:'2px 0 0 0'}}));
Map.add(legend);

// ── EXPORTS ────────────────────────────────────────────────────────────
Export.image.toDrive({image:contaminationOnly,description:'INFN8VZN_Mercury_Disturbance_2026',
  folder:'GalaSat',fileNamePrefix:'infn8vzn_mercury_disturbance_2026',
  region:AOI,scale:10,crs:'EPSG:4326',maxPixels:1e13,fileFormat:'GeoTIFF'});
Export.image.toDrive({image:turbidWater2025,description:'INFN8VZN_Turbid_Water_2025',
  folder:'GalaSat',fileNamePrefix:'infn8vzn_turbid_water_2025',
  region:AOI,scale:10,crs:'EPSG:4326',maxPixels:1e13,fileFormat:'GeoTIFF'});
Export.image.toDrive({image:flood_zone,description:'INFN8VZN_Flood_Zone',
  folder:'GalaSat',fileNamePrefix:'infn8vzn_flood_zone',
  region:AOI,scale:30,crs:'EPSG:4326',maxPixels:1e13,fileFormat:'GeoTIFF'});
print('── Exports submitted — check Tasks tab ──');
