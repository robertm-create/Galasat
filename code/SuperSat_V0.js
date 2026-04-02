// ============================================================
// SuperSat_V0 (GEE filename: SuperSat_V0 — submitted as "Galastat")
// MEMORY TAG: SuperSat_V0  | INTERNAL VERSION: prototype / "GalaSat baseline"
// SOURCE: claude.ai project "Technical - Repository"
//         chat "SuperSat_V1 code tagging" (Apr 20, 2026)
//         pasted attachment, 168 lines, 6.25 KB
// EXTRACTED: 2026-04-25
//
// LINEAGE: Older than V1/V2, newer than GalaSat_v2/v3.
//   Adds vs Gen1 baseline: GHSL urban mask, applyMinArea minimum-pixel filter, L8 detection layers.
//   Still missing: RF classifier, polygon assets, GMW proximity, turbid water, JRC flood,
//                  road masking, river/community layers.
//
// KNOWN ISSUE (carried into V1 + V2):
//   The mineS2 / mineL8 BSI denominator is written as
//     b11.add(b4).add(b8a.add(b2))
//   This adds (b8a + b2) as a sub-expression first, then sums with b11.add(b4).
//   Arithmetic-equivalent to (B11 + B4 + B8A + B2) in GEE's API but reads ambiguously.
//   Standardize as four independent additions in V3+ for clarity.
//
// NOTE: Whitespace partially collapsed during extraction. Logic and semicolons intact — runnable.
// ============================================================

var box = ee.Geometry.Rectangle([-2.3670, 5.8599, -1.4829, 6.5807]);
Map.centerObject(box, 10);
Map.addLayer(ee.Image().paint(box, 0, 2), {palette: ['FF0000']}, 'Boundary');

var ghsl = ee.ImageCollection('JRC/GHSL/P2016/SMOD_POP_GLOBE_V1')
  .filterBounds(box).mosaic().clip(box);
var urbanMask = ghsl.select('smod_code').lt(3);

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

function compL8(d1, d2) {
  return ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
    .filterBounds(box).filterDate(d1, d2)
    .map(maskL8).median().clip(box);
}

function compL8rgb(d1, d2) {
  return compL8(d1, d2)
    .select(['SR_B4','SR_B3','SR_B2'], ['B4','B3','B2']);
}

function applyMinArea(mask, minPixels) {
  var connected = mask.connectedPixelCount(200, true);
  return mask.updateMask(connected.gte(minPixels));
}

function mineS2(img) {
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
  var combined = a.or(b).or(c).or(d);
  var filtered = applyMinArea(combined, 12);
  var finalMask = filtered.and(valid).and(urbanMask);
  return bsi.updateMask(finalMask);
}

function mineL8(img) {
  var b2 = img.select('SR_B2');
  var b3 = img.select('SR_B3');
  var b4 = img.select('SR_B4');
  var b5 = img.select('SR_B5');
  var b6 = img.select('SR_B6');
  var valid = b2.mask().and(b3.mask()).and(b4.mask())
    .and(b5.mask()).and(b6.mask());
  var bsi = b6.add(b4).subtract(b5.add(b2)).divide(b6.add(b4).add(b5.add(b2)));
  var ndvi  = b5.subtract(b4).divide(b5.add(b4));
  var ndwi  = b3.subtract(b5).divide(b3.add(b5));
  var ndbi  = b6.subtract(b5).divide(b6.add(b5));
  var mndwi = b3.subtract(b6).divide(b3.add(b6));
  var a = bsi.gt(0.005).and(ndvi.lt(0.35)).and(ndbi.lt(0.05));
  var b = ndwi.gt(0.02).and(ndvi.lt(0.15)).and(bsi.gt(-0.08));
  var c = b4.gt(0.07).and(b3.gt(0.05)).and(ndvi.lt(0.1)).and(ndbi.lt(0.0));
  var d = mndwi.gt(0.0).and(ndvi.lt(0.1)).and(b4.gt(0.05));
  var combined = a.or(b).or(c).or(d);
  var filtered = applyMinArea(combined, 6);
  var finalMask = filtered.and(valid).and(urbanMask);
  return bsi.updateMask(finalMask);
}

function countHa(img, label, scale) {
  var area = img.mask().multiply(ee.Image.pixelArea());
  var total = area.reduceRegion({reducer: ee.Reducer.sum(), geometry: box, scale: scale, maxPixels: 1e10});
  var bandName = total.keys().get(0);
  print(label, ee.Number(total.get(bandName)).divide(10000), 'ha');
}

var tc   = {bands: ['B4', 'B3', 'B2'], min: 0.0, max: 0.3, gamma: 1.4};
var tcL8 = {bands: ['SR_B4', 'SR_B3', 'SR_B2'], min: 0.0, max: 0.3, gamma: 1.4};
var mv   = {min: 0.0, max: 0.3, palette: ['ffff00', 'fd8d3c', 'f03b20', 'bd0026']};

var l8_14 = compL8('2013-10-01', '2014-04-30');
var l8_15 = compL8('2014-10-01', '2015-04-30');
var s2_16 = compS2('2015-06-01', '2016-12-31');
var s2_18 = compS2('2017-06-01', '2018-06-30');
var s2_20 = compS2('2019-10-01', '2020-04-30');
var s2_22 = compS2('2021-10-01', '2022-04-30');
var s2_24 = compS2('2023-10-01', '2024-04-30');
var s2_25 = compS2('2024-10-01', '2025-04-30');
var s2_26 = compS2('2025-10-01', '2026-03-09');

var l8_16 = compL8rgb('2015-06-01', '2016-12-31');
var l8_18 = compL8rgb('2017-06-01', '2018-06-30');

var d16 = s2_16.select(['B4','B3','B2']).unmask(l8_16).clip(box);
var d18 = s2_18.select(['B4','B3','B2']).unmask(l8_18).clip(box);

Map.addLayer(l8_14, tcL8, '2014 True Color L8');
Map.addLayer(l8_15, tcL8, '2015 True Color L8');
Map.addLayer(d16,   tc,   '2016 True Color');
Map.addLayer(d18,   tc,   '2018 True Color');
Map.addLayer(s2_20, tc,   '2020 True Color');
Map.addLayer(s2_22, tc,   '2022 True Color');
Map.addLayer(s2_24, tc,   '2024 True Color');
Map.addLayer(s2_25, tc,   '2025 True Color');
Map.addLayer(s2_26, tc,   '2026 True Color');

var ml14 = mineL8(l8_14);
var ml15 = mineL8(l8_15);
var m16  = mineS2(s2_16);
var m18  = mineS2(s2_18);
var m20  = mineS2(s2_20);
var m22  = mineS2(s2_22);
var m24  = mineS2(s2_24);
var m25  = mineS2(s2_25);
var m26  = mineS2(s2_26);

Map.addLayer(ml14, mv, '2014 Mining L8');
Map.addLayer(ml15, mv, '2015 Mining L8');
Map.addLayer(m16,  mv, '2016 Mining');
Map.addLayer(m18,  mv, '2018 Mining');
Map.addLayer(m20,  mv, '2020 Mining');
Map.addLayer(m22,  mv, '2022 Mining');
Map.addLayer(m24,  mv, '2024 Mining');
Map.addLayer(m25,  mv, '2025 Mining');
Map.addLayer(m26,  mv, '2026 Mining');

print('GHSL urban mask applied');
print('Min area filter 12px S2 / 6px L8');
print('2014-2015 Landsat 8 30m');
print('2016-2026 Sentinel-2 10-20m');

countHa(ml14, '2014 Mining Ha L8', 30);
countHa(ml15, '2015 Mining Ha L8', 30);
countHa(m16,  '2016 Mining Ha S2', 20);
countHa(m18,  '2018 Mining Ha S2', 20);
countHa(m20,  '2020 Mining Ha S2', 20);
countHa(m22,  '2022 Mining Ha S2', 20);
countHa(m24,  '2024 Mining Ha S2', 20);
countHa(m25,  '2025 Mining Ha S2', 20);
countHa(m26,  '2026 Mining Ha S2', 20);

print('Done');
