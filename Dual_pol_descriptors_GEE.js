/*-----------------------------------------------------------------------------------------------------

This GEE script calculates dual-polarimetric descriptors 
(co-pol purity parameter mc, Pseudo scattering entropy Hc and psuedo scattering type parameter Theta_c) 
for Sentinel-1 GRD data.

OUTPUT: 

1) Temporal Sentinel-1 dialy mosaic scenes over the given ROI with following paramters as layers : 
   Hc, Theta_c, mc, class, ratio, VV, VH, inc 

2) Extracted values of the above paramters for given list of sampling points in *.csv format


Author Details: 
Narayana Rao B.
206-MRSLab, CSRE,
IIT Bombay, India.
email: bnarayanarao@iitb.ac.in
web: https://narayana-rao.github.io

A detailed explanation of the implemented algorithm can be found in the following articles.

Narayanarao Bhogapurapu, Subhadip Dey, Avik Bhattacharya, Dipankar Mandal, 
Juan Lopez-Sanchez, Heather McNairn, Carlos Lopez-Martinez and Y. S. Rao 2021
“Dual-polarimetric descriptors from Sentinel-1 GRD SAR data for crop growth assessment”. 
ISPRS Journal of Photogrammetry and Remote Sensing. 20-35, 178. 
doi: 10.1016/j.isprsjprs.2021.05.013

Narayanarao Bhogapurapu, Subhadip Dey, Dipankar Mandal, Avik Bhattacharya, 
L. Karthikeyan, Heather McNairn and Y. S. Rao 2022
“Soil Moisture Retrieval Over Croplands Using dual-pol L-band GRD SAR Data”. 
Remote Sensing of Environment. Volume 271, 2022, Pages 112900, ISSN 0034-4257 
doi: 10.1016/j.rse.2022.112900
-------------------------------------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------------------------

                      1) Import ROI(extent) and sampling points

----------------------------------------------------------------------------------------------*/


var extent = 
    ee.Geometry.Polygon(
        [[[-97.33036743800388, 49.68259937829183],
          [-97.33036743800388, 49.55804339464324],
          [-97.20951782862888, 49.55804339464324],
          [-97.20951782862888, 49.68259937829183]]], null, false);
          
var sample_pts = ee.FeatureCollection([
ee.Feature(ee.Geometry.Point(-98.04639258,49.68454278), {label: 'P1'}),
ee.Feature(ee.Geometry.Point(-98.04642806,49.68251043), {label: 'P2'}),
ee.Feature(ee.Geometry.Point(-98.04643989,49.68183298), {label: 'P3'}),
ee.Feature(ee.Geometry.Point(-98.04646585,49.68116492), {label: 'P4'}),

]);


/*----------------------------------------------------------------------------------------------

                      2) Cloud filtering and data preparation

----------------------------------------------------------------------------------------------*/


var ref_start=ee.Date('2016-08-15');
var ref_end = ee.Date('2016-09-30');

var window_size = 2.5; //window size for filtering
print('window size',window_size*2);

var S1 = ee.ImageCollection('COPERNICUS/S1_GRD')
        .filterDate(ref_start, ref_end)
        .filter(ee.Filter.eq('instrumentMode', 'IW'))
        .filter(ee.Filter.eq('orbitProperties_pass', 'ASCENDING')) 
        .select('VV','VH','angle')
        .sort('system:time_start', false)
        .filterBounds(extent);


// Difference in days between start and finish
var diff = ref_end.difference(ref_start, 'day');

// Make a list of all dates
var range = ee.List.sequence(0, diff.subtract(1)).map(function(day)
{return ref_start.advance(day,'day')})

// Funtion for iteraton over the range of dates
var day_mosaics = function(date, newlist) {
  // Cast
  date = ee.Date(date)
  newlist = ee.List(newlist)

  // Filter collection between date and the next day
  var filtered = S1.filterDate(date, date.advance(1,'day'))

  // Make the mosaic
  var image = ee.Image(filtered.mosaic());
  // copy image meta
  image = image
          .set('system:time_start', filtered.first().get('system:time_start'))
          .set('system:index', filtered.first().get('system:index'))
          .set('system:id', filtered.first().get('system:id'))
          .set('system:version', filtered.first().get('system:version'))
          // .set('system:bands', filtered.first().get('system:bands'))
          .set('system:footprint', filtered.first().get('system:footprint'))
          ;

  // Add the mosaic to a list only if the collection has images
  return ee.List(ee.Algorithms.If(filtered.size(), newlist.add(image), newlist))
}

// Iterate over the range to make a new list, and then cast the list to an imagecollection
var newS1col = ee.ImageCollection(ee.List(range.iterate(day_mosaics, ee.List([]))))


/*----------------------------------------------------------------------------------------------

                      3) Generating Dual-pol descriptors and the clusters

----------------------------------------------------------------------------------------------*/


var m = newS1col.map(function(image) {
    var C11_mean = image.expression( '10 ** (VV / 10)', {'VV': image.select('VV')})
                  .reduceNeighborhood({
                    reducer: ee.Reducer.mean(),
                    kernel: ee.Kernel.square(window_size)
                    });
    var C22_mean = image.expression( '10 ** (VH / 10)', {'VH': image.select('VH')})
                  .reduceNeighborhood({
                    reducer: ee.Reducer.mean(),
                    kernel: ee.Kernel.square(window_size)
                    });

    var span = C11_mean.add(C22_mean);
    var ratio = C22_mean.divide(C11_mean);
    var vmask = C11_mean.subtract(C22_mean);
    vmask = vmask.expression('b(0) >0? 1:0');
    
    var m = (C11_mean.subtract(C22_mean).abs()).divide(span);
    var d_dpol = m.multiply(m).subtract(1).multiply(-1);
    var theta_c = ((C11_mean.subtract(C22_mean).abs()).multiply(span).multiply(m))
                        .divide((C11_mean.multiply(C22_mean)).add(span.pow(2).multiply(m.pow(2))))
                        .atan();
    theta_c = theta_c.multiply(180).divide(Math.PI);

    var p1 = C11_mean.divide(span);
    var p2 = C22_mean.divide(span);
    var cnst = ee.Number(2);
    var Hp1 = p1.multiply(p1.log10()).divide(cnst.log10()).multiply(-1);
    var Hp2 = p2.multiply(p2.log10()).divide(cnst.log10()).multiply(-1);
    var H = Hp1.add(Hp2);
    var q = ratio;
    var DpRVIc_n = q.multiply(q.add(ee.Number(3)));
    var DpRVIc_d = (q.add(ee.Number(1))).multiply(q.add(ee.Number(1)));
    var DpRVIc = DpRVIc_n.divide(DpRVIc_d);


    var H_rc = H.expression('b(0) >0 && b(0)<0.3 ? 1 : b(0) > 0.3 && b(0) <0.5? 2 : b(0)>0.5&&b(0) < 0.7  ? 3 :b(0)>0.7 && b(0)<1.0 ? 4: 0');
    var theta_c_rc = theta_c.expression('b(0)>0.0 && b(0) <15 ? 5 :  b(0)>15 &&  b(0)<30 ? 6 : b(0)>30 && b(0) < 45? 7 : 0');
    var C11_mean_db = C11_mean.log10().multiply(10);//Linear to dB conversion
    var C11_rc = C11_mean_db.expression('b(0)<-17?0:1'); // masking low dB returns (water)
    
    var out = H_rc.multiply(theta_c_rc).multiply(C11_rc);
    var out_rc = out.expression('b(0) ==7 ? 1 : b(0) == 14 ? 2 : b(0) == 21 ? 3: b(0) == 20 ? 6: b(0) == 24 ? 5: b(0) == 28 ? 4 : 0');
    
    
    //Masked values
    m = (m.updateMask(vmask)).updateMask(C11_rc);
    H=(H.updateMask(vmask)).updateMask(C11_rc);
    theta_c=(theta_c.updateMask(vmask)).updateMask(C11_rc);
    DpRVIc=(DpRVIc.updateMask(vmask)).updateMask(C11_rc);
    out_rc=(out_rc.updateMask(vmask)).updateMask(C11_rc);
    ratio=(ratio.updateMask(vmask)).updateMask(C11_rc);

    var out_raster = H.addBands([theta_c.select('constant_mean'),
                                m.select('constant_mean'),
                                
                                out_rc.select('constant').toDouble(),
                                ratio.select('constant_mean'),
                                C11_mean.select('constant_mean'),
                                C22_mean.select('constant_mean'),
                                DpRVIc.select('constant_mean'),
                                image.select('angle')]);

    out_raster = out_raster.select(
        ['constant_mean', 'constant_mean_1','constant_mean_1_1','constant','constant_mean_2','constant_mean_3','constant_mean_4','constant_mean_5','angle'], // old names
        ['Hc', 'Theta_c','mc','class','ratio','VV','VH','DpRVIc','inc']             
        );
    return out_raster.set('system:time_start', image.get('system:time_start'));

});


// output visualization

var jet_cmap = [' #000080 ', ' #0000bd ', ' #0000fa ', ' #0022ff ', ' #0057ff ', ' #008dff ',
' #00c3ff ', ' #0ff8e8 ', ' #3affbc ', ' #66ff91 ', ' #91ff66 ', ' #bcff3a ', ' #e8ff0f ', ' #ffd500 ',
' #ffa400 ', ' #ff7200 ', ' #ff4000 ', ' #fa0e00 ', ' #bd0000 ', ' #800000 ',]

Map.centerObject(extent,15);
Map.addLayer(ee.Image(m.select('Hc').first()),{min:0,max:1,palette:jet_cmap},'Hc');
Map.addLayer(ee.Image(m.select('mc').first()),{min:0,max:1,palette:jet_cmap},'mc');
Map.addLayer(ee.Image(m.select('Theta_c').first()),{min:0,max:45,palette:jet_cmap},'Theta_c');
Map.addLayer(ee.Image(m.select('DpRVIc').first()),{min:0,max:1,palette:jet_cmap},'DpRVIc');

/*----------------------------------------------------------------------------------------------

                      4) Exporting the data in ratser format and csv

----------------------------------------------------------------------------------------------*/


var bandcol = ee.List(['Hc','Theta_c','mc','DpRVIc','class','ratio','VV','VH','inc']); 
var bandsize = bandcol.size().getInfo();
for (var i = 0; i < bandsize; i++) {
var band = ee.String(bandcol.get(i));
var sample_pts = sample_pts.map(function(feature) {
  return ee.Feature(feature.geometry(), {'id': feature.id()})
});

var triplets = m.map(function(image) {
  return image.select(band).reduceRegions({
    collection: sample_pts, 
    reducer: ee.Reducer.first().setOutputs([band]), 
    scale: 30,
  }).map(function(feature) {
    var dpgrd = ee.List([feature.get(band), -9999])
      .reduce(ee.Reducer.firstNonNull())
    return feature.set({band : dpgrd, 'imageID': image.id()})
    })
  }).flatten();

var format = function(table, rowId, colId) {
  var rows = table.distinct(rowId); 
  var joined = ee.Join.saveAll('matches').apply({
    primary: rows, 
    secondary: table, 
    condition: ee.Filter.equals({
      leftField: rowId, 
      rightField: rowId
    })
  });
        
  return joined.map(function(row) {
      var values = ee.List(row.get('matches'))
        .map(function(feature) {
          feature = ee.Feature(feature);
          return [feature.get(colId), feature.get(band)];
        });
      return row.select([rowId]).set(ee.Dictionary(values.flatten()));
    });
};

var sentinelResults = format(triplets, 'id', 'imageID');

var merge = function(table, rowId) {
  return table.map(function(feature) {
    var id = feature.get(rowId)
    var allKeys = feature.toDictionary().keys().remove(rowId)
    var substrKeys = ee.List(allKeys.map(function(val) { 
        return ee.String(val).slice(0,8)}
        ))
    var uniqueKeys = substrKeys.distinct()
    var pairs = uniqueKeys.map(function(key) {
      var matches = feature.toDictionary().select(allKeys.filter(ee.Filter.stringContains('item', key))).values()
      var val = matches.reduce(ee.Reducer.max())
      return [key, val]
    })
    return feature.select([rowId]).set(ee.Dictionary(pairs.flatten()))
  })
}
var sentinelMerged = merge(sentinelResults, 'id');
// print(ee.String(band));
var band = bandcol.get(i);
Export.table.toDrive({
    collection: sentinelResults,
    description: bandcol.get(i)+'_time_series',
    folder: 'dpgrd_out',
    fileNamePrefix:band+'_time_series',
    fileFormat: 'CSV'
});
}


var ExportCol = function(col, folder, scale, type,
                         nimg, maxPixels, region) {
    type = type || "float";
    nimg = nimg || 500;
    scale = scale || 30;
    maxPixels = maxPixels || 1e12;

    var colList = col.toList(nimg);
    var n = colList.size().getInfo();

    for (var i = 0; i < n; i++) {
      var img = ee.Image(colList.get(i));
      var id = img.id().getInfo();
      region = region || img.geometry().bounds().getInfo()["coordinates"];

      var imgtype = {"float":img.toFloat(), 
                     "byte":img.toByte(), 
                     "int":img.toInt(),
                     "double":img.toDouble()
                    }

      Export.image.toDrive({
        image:imgtype[type],
        description: id,
        folder: folder,
        fileNamePrefix: id,
        region: region,
        scale: scale,
        maxPixels: maxPixels})
    }
  }

//Uncomment the below line to export all the avaialble scenes and corresponding descriptors in the Geotiff format. 

// ExportCol(m, 'dpgrd_out', 30,"double",100,1e12,extent) 


/*----------------------------------------------------------------------------------------------

                                    END of the script

----------------------------------------------------------------------------------------------*/
