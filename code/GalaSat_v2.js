// GalaSat v2 — Last Confirmed Working Build
// INFN8VZN Fed Inc. | Federal Corporation No. 1750155-2
// Study Area: Ashanti Region, Ghana — Highest concentrated galamsey district
// Seven time periods: 2016, 2018, 2020, 2022, 2024, 2025, 2026
// Sentinel-2 primary, Landsat 8 fills gaps for 2016 and 2018 display only
// Mining detection: pure Sentinel-2 with valid pixel mask
// Hectare calculation: printed to console per year
// Resolution: Sentinel-2 10-20m
// Last confirmed working: 2026-03-09
//
// MEMORY TAG: GalaSat_v2  | INTERNAL VERSION: v2 (Gen 1 baseline, canonical reference)
// SOURCE: claude.ai project "Technical - Repository"
//         chat "SuperSat_V1 code tagging" (Apr 20, 2026)
//         pasted attachment, 140 lines
// EXTRACTED: 2026-04-25
//
// CONFIRMED CONSOLE OUTPUT:
// 2016 Mining Hectares: 39879.70 ha (NOTE: may include artifacts — see methodology)
// 2018 Mining Hectares: 7418.15 ha (Operation Vanguard suppression period)
// 2020 Mining Hectares: 17929.27 ha
// 2022 Mining Hectares: 23606.24 ha
// 2024 Mining Hectares: 30686.04 ha
// 2025 Mining Hectares: 39824.11 ha
// 2026 Mining Hectares: 43676.67 ha
//
// LINEAGE NOTE: This is the unfiltered ceiling. No GHSL urban mask, no applyMinArea, no L8 detection.
// Every later build adds suppression — should produce LOWER hectare numbers. If a future build
// approaches these figures, something broke upstream. These are the regression canary.
//
// NOTE: Whitespace partially collapsed during extraction. Logic and semicolons intact — runnable.

var box = ee.Geometry.Rectangle([-2.3670, 5.8599, -1.4829, 6.5807]);
Map.centerObject(box, 10);
Map.addLayer(ee.Image().paint(box, 0, 2), {palette: ['FF0000']}, 'Boundary');

function maskS2(img) {
  var qa = img.select('QA60');
  var qm = qa.bitwiseAnd(1 << 10).eq(0).and(qa.bitwiseAnd(1 << 11).eq(0));
  var scl = img.select('SCL');
  var sm = scl.neq(3).and(scl.neq(8)).and(scl.neq(9)).and(scl.neq(10)).and(scl.neq(11));
  return img.updateMask(qm).updateMask(sm).divide(10000);
}

function maskL8(img) {
  var qa = img.select('QA_PIXEL');
  var m = qa.bitwiseAnd(1 << 3).eq(0).and(qa.bitwiseAnd(1 << 4).eq(0));
  return img.updateMask(m).multiply(0.0000275).add(-0.2);
}

function compS2(d1, d2) {
  return ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
    .filterBounds(box).filterDate(d1, d2)
    .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 30))
    .map(maskS2).median().clip(box);
}

function compL8rgb(d1, d2) {
  return ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
    .filterBounds(box).filterDate(d1, d2)
    .map(maskL8).median()
    .select(['SR_B4','SR_B3','SR_B2'], ['B4','B3','B2'])
    .clip(box);
}

function mine(img) {
  var b2 = img.select('B2');
  var b3 = img.select('B3');
  var b4 = img.select('B4');
  var b8 = img.select('B8');
  var b8a = img.select('B8A');
  var b11 = img.select('B11');
  var valid = b2.mask().and(b3.mask()).and(b4.mask())
    .and(b8.mask()).and(b8a.mask()).and(b11.mask());
  var bsi = b11.add(b4).subtract(b8a.add(b2)).divide(b11.add(b4).add(b8a.add(b2)));
  var ndvi  = b8.subtract(b4).divide(b8.add(b4));
  var ndwi  = b3.subtract(b8).divide(b3.add(b8));
  var ndbi  = b11.subtract(b8).divide(b11.add(b8));
  var mndwi = b3.subtract(b11).divide(b3.add(b11));
  var a = bsi.gt(0.005).and(ndvi.lt(0.35)).and(ndbi.lt(0.05));
  var b = ndwi.gt(0.02).and(ndvi.lt(0.15)).and(bsi.gt(-0.08));
  var c = b4.gt(0.07).and(b3.gt(0.05)).and(ndvi.lt(0.1)).and(ndbi.lt(0.0));
  var d = mndwi.gt(0.0).and(ndvi.lt(0.1)).and(b4.gt(0.05));
  return bsi.updateMask(a.or(b).or(c).or(d)).updateMask(valid);
}

function countHectares(img, label) {
  var binary = img.mask().rename('mining');
  var area = binary.multiply(ee.Image.pixelArea()).rename('area');
  var total = area.reduceRegion({reducer: ee.Reducer.sum(), geometry: box, scale: 20, maxPixels: 1e10});
  var ha = ee.Number(total.get('area')).divide(10000);
  print(label, ha, 'hectares');
}

var tc = {bands: ['B4', 'B3', 'B2'], min: 0.0, max: 0.3, gamma: 1.4};
var mv = {min: 0.0, max: 0.3, palette: ['ffff00', 'fd8d3c', 'f03b20', 'bd0026']};

// S2 composites
var s2_16 = compS2('2015-06-01', '2016-12-31');
var s2_18 = compS2('2017-06-01', '2018-06-30');
var s2_20 = compS2('2019-10-01', '2020-04-30');
var s2_22 = compS2('2021-10-01', '2022-04-30');
var s2_24 = compS2('2023-10-01', '2024-04-30');
var s2_25 = compS2('2024-10-01', '2025-04-30');
var s2_26 = compS2('2025-10-01', '2026-03-09');  // v2 marker (v3 has 2026-03-08)

// L8 fills for 2016 and 2018 display only
var l8_16 = compL8rgb('2015-06-01', '2016-12-31');
var l8_18 = compL8rgb('2017-06-01', '2018-06-30');

// Display layers
var d16 = s2_16.select(['B4','B3','B2']).unmask(l8_16).clip(box);
var d18 = s2_18.select(['B4','B3','B2']).unmask(l8_18).clip(box);

Map.addLayer(d16, tc, '2016 True Color');
Map.addLayer(d18, tc, '2018 True Color');
Map.addLayer(s2_20, tc, '2020 True Color');
Map.addLayer(s2_22, tc, '2022 True Color');
Map.addLayer(s2_24, tc, '2024 True Color');
Map.addLayer(s2_25, tc, '2025 True Color');
Map.addLayer(s2_26, tc, '2026 True Color');

// Mining detection — pure S2 only
var m16 = mine(s2_16);
var m18 = mine(s2_18);
var m20 = mine(s2_20);
var m22 = mine(s2_22);
var m24 = mine(s2_24);
var m25 = mine(s2_25);
var m26 = mine(s2_26);

Map.addLayer(m16, mv, '2016 Mining');
Map.addLayer(m18, mv, '2018 Mining');
Map.addLayer(m20, mv, '2020 Mining');
Map.addLayer(m22, mv, '2022 Mining');
Map.addLayer(m24, mv, '2024 Mining');
Map.addLayer(m25, mv, '2025 Mining');
Map.addLayer(m26, mv, '2026 Mining');

countHectares(m16, '2016 Mining Hectares');
countHectares(m18, '2018 Mining Hectares');
countHectares(m20, '2020 Mining Hectares');
countHectares(m22, '2022 Mining Hectares');
countHectares(m24, '2024 Mining Hectares');
countHectares(m25, '2025 Mining Hectares');
countHectares(m26, '2026 Mining Hectares');

print('Done');
