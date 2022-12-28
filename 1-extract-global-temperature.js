// The javascript code can be accessed through Google Earth Engine cloud computing platform at:
// https://code.earthengine.google.com/3596a7d8f40f4195cf8102a0ce81ac3c

// The code is as follows:

var myd11a1 = ee.ImageCollection("MODIS/006/MYD11A1"),
    myd14a1 = ee.ImageCollection("MODIS/006/MYD14A1"),
    myd13a2 = ee.ImageCollection("MODIS/006/MYD13A2"),
    mcd12q1 = ee.ImageCollection("MODIS/006/MCD12Q1");

var roi = /* color: #bf04c2 */ee.Geometry.Polygon(
        [[[-180, 90],
          [0, 90],
          [180, 90],
          [180, -90],
          [0, -90],
          [-180, -90]]], null, false);

var ndviMask = myd13a2.filterDate("2015-1-1", "2016-1-1")
                      .select("NDVI")
                      .mosaic()
                      .mask();

// var roi = /* color: #0b4a8b */ee.Geometry.Point([-77.32248808,	-9.320021073]);

// extract daily maximym and max/min lst for earch location
// mid&high/high fire excluded
// QA is good/not good
// After mosaicing

// Funtion for iteraton over the range of dates
function day_mosaics(date, newlist, rawImgCol) {
  // Cast
  date = ee.Date(date);
  newlist = ee.List(newlist);
  // Filter collection between date and the next day
  var filtered = rawImgCol.filterDate(date, date.advance(1,'day'));
  // Make the mosaic
  var image = ee.Image(filtered.mosaic());
  image = image.set("system:time_start", date.format("yyyyMMdd"));
  image = image.set("system:index", date.format("yyyy_MM_dd"));
  // Add the mosaic to a list only if the collection has images
  return ee.List(ee.Algorithms.If(filtered.size(), newlist.add(image), newlist));
}

function getGoodLSTByBand(startDateStr, endDateStr, qTag, bandTag) {
  var startDate = ee.Date(startDateStr);
  var endDate = ee.Date(endDateStr);
  var days = endDate.difference(startDate, "day");
  var dayList = ee.List.sequence(0, days.subtract(1));
  var qc_band_name = "QC_"+bandTag;
  var lst_band_name = "LST_"+bandTag+"_1km";
  var imgList = dayList.map(function(day){
    day = ee.Number(day);
    var curDate = startDate.advance(day, "day");
    var nextDate = curDate.advance(1, "day");
    var img = myd11a1.filterDate(curDate, nextDate)
                    .select(["LST_"+bandTag+"_1km", "QC_"+bandTag])
                    .first();
                   
    //only best QA value/all values without quality control
    var lstImg = null;
    var constant_img = ee.Image.constant([0, 0])
                         .select(["constant_0","constant_1"], [lst_band_name, qc_band_name])
                         .toByte();
    constant_img = constant_img.updateMask(constant_img.select(qc_band_name));
    if (qTag === 0) {
      lstImg = ee.Algorithms.If(
        img,
        img,
        constant_img
      );
    } else {
      
      lstImg = ee.Algorithms.If(
        // Good Quality
        img,
        // img.updateMask(img.select(qc_band_name).eq(0)),
        img.updateMask(img.select(qc_band_name).bitwiseAnd(1<<1).eq(0)
                          .and(img.select(qc_band_name).bitwiseAnd(1<<2).eq(0))
                          .and(img.select(qc_band_name).bitwiseAnd(1<<3).eq(0))
                          .and(img.select(qc_band_name).lt(1<<5))),
        constant_img
      );
    }
    ///////////end/////////////
    
    lstImg = ee.Image(lstImg);
    lstImg = lstImg.set("system:time_start", curDate.millis());
    lstImg = lstImg.set("system:index", curDate.format("yyyy_MM_dd"));
    return lstImg;
  });
  var imgCol = ee.ImageCollection(imgList);
  return imgCol;
}

function getGoodLST(startDateStr, endDateStr, qTag) {
  var startDate = ee.Date(startDateStr);
  var endDate = ee.Date(endDateStr);
  var days = endDate.difference(startDate, "day");
  var dayList = ee.List.sequence(0, days.subtract(1));
  var imgList = dayList.map(function(day){
    day = ee.Number(day);
    var curDate = startDate.advance(day, "day");
    var nextDate = curDate.advance(1, "day");
    var img = myd11a1.filterDate(curDate, nextDate)
                    .select(["LST_Day_1km", "QC_Day", "LST_Night_1km", "QC_Night"])
                    .first();
                   
    //only best QA value/all values without quality control
    var lstImg = null;
    var constant_img = ee.Image.constant([0, 0, 0, 0])
                         .select(
                           ["constant_0","constant_1","constant_2","constant_3"], 
                           ["LST_Day_1km", "QC_Day", "LST_Night_1km", "QC_Night"])
                         .toByte();
    constant_img = constant_img.updateMask(constant_img.select("QC_Day"));
    if (qTag === 0) {
      
      lstImg = ee.Algorithms.If(
        img,
        img,
        constant_img
      );
    } else {
     
      lstImg = ee.Algorithms.If(
        // Good Quality
        img,
        // img.updateMask(img.select("QC_Day").eq(0).and(img.select("QC_Night").eq(0))),
        img.updateMask(img.select("QC_Day").bitwiseAnd(1<<1).eq(0)
                          .and(img.select("QC_Day").bitwiseAnd(1<<2).eq(0))
                          .and(img.select("QC_Day").bitwiseAnd(1<<3).eq(0))
                          .and(img.select("QC_Day").lt(1<<5))
                          .and(
                            img.select("QC_Night").bitwiseAnd(1<<1).eq(0)
                               .and(img.select("QC_Night").bitwiseAnd(1<<2).eq(0))
                               .and(img.select("QC_Night").bitwiseAnd(1<<3).eq(0))
                               .and(img.select("QC_Night").lt(1<<5))
                          )),
        constant_img
      );
    }
    ///////////end/////////////
    
    lstImg = ee.Image(lstImg);
    lstImg = lstImg.set("system:time_start", curDate.millis());
    lstImg = lstImg.set("system:index", curDate.format("yyyy_MM_dd"));
    return lstImg;
  });
  var imgCol = ee.ImageCollection(imgList);
  return imgCol;
}

function getMosaicLst(imgCol, startDate, endDate) {
  var lst = imgCol.filter(ee.Filter.date(startDate, endDate));
  var diff = endDate.difference(startDate, 'day');
  // Make a list of all dates
  var range = ee.List.sequence(0, diff.subtract(1))
                .map(function(day){
                  return startDate.advance(day,'day');
                }); 
  // Iterate over the range to make a new list, and then cast the list to an imagecollection
  var lstList = range.iterate(function(date, newlist) {
    return day_mosaics(date, newlist, lst);
  }, ee.List([]));
  var MosaicLst = ee.ImageCollection(ee.List(lstList));
  return MosaicLst;
}

function getMosaicFire(startDate, endDate) {
  var fire = myd14a1.filter(ee.Filter.date(startDate, endDate))
                    .select('FireMask');

  var diff = endDate.difference(startDate, 'day');
  // Make a list of all dates
  var range = ee.List.sequence(0, diff.subtract(1))
                .map(function(day){
                  return startDate.advance(day,'day');
                });
  var fireList = range.iterate(function(date, newlist) {
    return day_mosaics(date, newlist, fire);
  }, ee.List([]));
  // Iterate over the range to make a new list, and then cast the list to an imagecollection
  var MosaicFire = ee.ImageCollection(ee.List(fireList));
  return MosaicFire;
}

function getNDVIImage(start_date, end_date) {
  var sMillis = ee.Date(start_date).millis();
  var eMillis = ee.Date(end_date).millis();
  var ndviImage = myd13a2.filter(
                    ee.Filter.or(
                      ee.Filter.and(
                        ee.Filter.lte("system:time_start", sMillis),
                        ee.Filter.gte("system:time_end", sMillis)
                      ),
                      ee.Filter.and(
                        ee.Filter.gte("system:time_start", sMillis),
                        ee.Filter.lte("system:time_end", eMillis)
                      ),
                      ee.Filter.and(
                        ee.Filter.lte("system:time_start", eMillis),
                        ee.Filter.gte("system:time_end", eMillis)
                      )
                    ))
                    .select("NDVI")
                    .map(function(image) {
                      return image.multiply(0.0001);
                    })
                    .first();
  ndviImage = ee.Algorithms.If(ndviImage, ndviImage, ee.Image.constant(0).updateMask(ndviMask).toByte().rename("NDVI"));
  ndviImage = ee.Image(ndviImage);
  return ndviImage;
}

function fireMask(image, MosaicFire) {
  var index = image.get("system:index");
  var fireImage = MosaicFire.filter(ee.Filter.eq("system:index", index))
                            .first();
  var lowFire = 7;                          
  var midFire = 8;
  var highFire = 9;
  var maskImage = ee.Image.constant(1).toByte();
  var mask = ee.Algorithms.If(
    fireImage, 
    fireImage.select("FireMask").eq(midFire)
            .or(fireImage.select("FireMask").eq(highFire))
            .or(fireImage.select("FireMask").eq(lowFire)), 
    // fireImage.select("FireMask").eq(highFire),      
    ee.Image.constant(0).toByte()
  );
  maskImage = maskImage.where(ee.Image(mask), 0);
  return image.updateMask(maskImage);
}

function getMosaicIGBP(land_covers, start_date, end_date) {
  var syear = ee.Date(start_date).get("year");
  var eyear = ee.Date(end_date).get("year");
  var igbp = mcd12q1.filterDate(ee.Date.fromYMD(syear, 1, 1), ee.Date.fromYMD(ee.Number(eyear).add(1), 1, 1))
                    .select("Land_Cover_Type_1")
                    .map(function(image) {
                      land_covers = ee.List(land_covers);
                      var maskImgList = land_covers.map(function(land_cover){
                        land_cover = ee.Number(land_cover);
                        var mask = image.eq(land_cover);
                        return mask;
                      });
                      var maskImgCol = ee.ImageCollection.fromImages(maskImgList);
                      var maskImg = maskImgCol.sum();
                      image = image.updateMask(maskImg);
                      var year = ee.Date(image.get("system:time_start")).get("year");
                      image = image.set("year", ee.Number(year).toInt());
                      return image;
                    });
  return igbp;
}

function igbpMask(image, MosaicIGBP) {
  var time_start = ee.Date.parse("yyyyMMdd", image.get("system:time_start"));
  var year = ee.Number(time_start.get("year")).toInt();
  var igbpImage = MosaicIGBP.filter(ee.Filter.eq("year", year))
                            .first();
  
  var mask = ee.Algorithms.If(
    igbpImage, 
    igbpImage,      
    ee.Image.constant(1).toByte()
  );
  mask = ee.Image(mask);
  return image.updateMask(mask);
}

function scaleDayLST(image) {
  var time_start = image.get("system:time_start");
  var qc_band = image.select("QC_Day");
  image = image.select("LST_Day_1km")
               .multiply(0.02)
               .subtract(273.15)
               .toDouble()
               .addBands(qc_band);
  image = image.set("system:time_start", time_start);
  return image; 
}

function scaleNightLST(image) {
  var time_start = image.get("system:time_start");
  var qc_band = image.select("QC_Night");
  image = image.select("LST_Night_1km")
               .multiply(0.02)
               .subtract(273.15)
               .toDouble()
               .addBands(qc_band);
  image = image.set("system:time_start", time_start);
  return image; 
}

function scaleLST(image) {
  var time_start = image.get("system:time_start");
  var qc_night = image.select("QC_Night");
  var qc_day = image.select("QC_Day");
  image = image.select(["LST_Day_1km", "LST_Night_1km"])
               .multiply(0.02)
               .subtract(273.15)
               .toDouble()
               .addBands(qc_day)
               .addBands(qc_night);
  image = image.set("system:time_start", time_start);
  return image; 
}

/**
 * 
 * */
function getDayLSTImages(start, finish, qTag, ndvi_threshold, MosaicFire) {
  var LST = getGoodLSTByBand(start, finish, qTag, "Day");
  var MosaicLst = getMosaicLst(LST, start, finish);
  var imgCol = null;
  if (ndvi_threshold === null) {
    imgCol = MosaicLst.map(function(image) {
                        return fireMask(image, MosaicFire);
                      })
                      .map(scaleDayLST);
  } else {
    imgCol = MosaicLst.map(function(image) {
                        var time_start = ee.Date.parse("yyyyMMdd", image.get("system:time_start"));
                        var ndviImage = getNDVIImage(time_start, time_start.advance(1, "day"));
                        return image.updateMask(ndviImage.lte(ndvi_threshold));
                      })
                      .map(function(image) {
                        return fireMask(image, MosaicFire);
                      })
                      .map(scaleDayLST);
  }
  return imgCol;
}

/**
 * 
 * */
function getNightLSTImages(start, finish, qTag) {
  var LST = getGoodLSTByBand(start, finish, qTag, "Night");
  var MosaicLst = getMosaicLst(LST, start, finish);
  var imgCol = MosaicLst.map(scaleNightLST);
  return imgCol;
}

/**
 * 
 * */
function getAllLSTImages(start, finish, qTag, ndvi_threshold, MosaicFire) {
  var LST = getGoodLST(start, finish, qTag);
  var MosaicLst = getMosaicLst(LST, start, finish);
  var imgCol = null;
  if (ndvi_threshold === null) {
    imgCol = MosaicLst.map(function(image) {
                        return fireMask(image, MosaicFire);
                      })
                      .map(scaleLST);
  } else {
    imgCol = MosaicLst.map(function(image) {
                        var time_start = ee.Date.parse("yyyyMMdd", image.get("system:time_start"));
                        var ndviImage = getNDVIImage(time_start, time_start.advance(1, "day"));
                        return image.updateMask(ndviImage.lte(ndvi_threshold));
                      })
                      .map(function(image) {
                        return fireMask(image, MosaicFire);
                      })
                      .map(scaleLST);
  }
  
  return imgCol;
}

//Cal Max LST
function calculateMaxLST(image) {
  image = image.select(["LST_Day_1km", "QC_Day"])
               .addBands(ee.Image.pixelLonLat());
  var dict = image.reduceRegion({
    reducer: ee.Reducer.max(4)
               .setOutputs(["max","max_qc","max_lon","max_lat"]),
    geometry: roi, 
    scale: 1000,
    maxPixels: 1e13
  });
  image = image.set("max", ee.Number(dict.get("max")));
  image = image.set("max_qc", ee.Number(dict.get("max_qc")));
  image = image.set("max_lon", ee.Number(dict.get("max_lon")));
  image = image.set("max_lat", ee.Number(dict.get("max_lat")));
  return image;
}

//Cal Min LST
function calculateMinLST(image) {
  image = image.select(["LST_Night_1km", "QC_Night"])
               .addBands(ee.Image.pixelLonLat());
  var dict = image.reduceRegion({
    reducer: ee.Reducer.min(4)
               .setOutputs(["min","min_qc","min_lon","min_lat"]),
    geometry: roi, 
    scale: 1000,
    maxPixels: 1e13
  });
  image = image.set("min", ee.Number(dict.get("min")));
  image = image.set("min_qc", ee.Number(dict.get("min_qc")));
  image = image.set("min_lon", ee.Number(dict.get("min_lon")));
  image = image.set("min_lat", ee.Number(dict.get("min_lat")));
  return image;
}

//Cal LST Dif
function calculateMaxDifferLST(image) {
  var qc_day = image.select("QC_Day");
  var qc_night = image.select("QC_Night");
  var lst_day = image.select("LST_Day_1km");
  var lst_night = image.select("LST_Night_1km");
  var differImage = lst_day.subtract(lst_night)
                          .abs()
                          .rename("differ");
  differImage = differImage.addBands(lst_day)
                          .addBands(lst_night)
                          .addBands(qc_day)
                          .addBands(qc_night)
                          .addBands(ee.Image.pixelLonLat());
  var dict = differImage.reduceRegion({
    reducer: ee.Reducer.max(7)
               .setOutputs(["differ","lst_day", "lst_night","qc_day","qc_night","differ_lon","differ_lat"]),
    geometry: roi, 
    scale: 1000,
    maxPixels: 1e13
  });
  image = image.set("differ", ee.Number(dict.get("differ")));
  image = image.set("lst_day", ee.Number(dict.get("lst_day")));
  image = image.set("lst_night", ee.Number(dict.get("lst_night")));
  image = image.set("qc_day", ee.Number(dict.get("qc_day")));
  image = image.set("qc_night", ee.Number(dict.get("qc_night")));
  image = image.set("differ_lon", ee.Number(dict.get("differ_lon")));
  image = image.set("differ_lat", ee.Number(dict.get("differ_lat")));
  return image;
}


function getDailyMaxLST(start, finish, qTag, ndvi_threshold, MosaicFire, exportName,foldername) {
  var imgCol = getDayLSTImages(start, finish, qTag, ndvi_threshold, MosaicFire);
  imgCol = imgCol.map(calculateMaxLST);  
  Export.table.toDrive({
    collection: imgCol, 
    description: "Drive-"+exportName,
    fileNamePrefix: exportName,
    folder: foldername,
    fileFormat: "CSV"
  });
  return imgCol;
}


function getDailyMinLST(start, finish, qTag, exportName,foldername) {
  var imgCol = getNightLSTImages(start, finish, qTag);
  imgCol = imgCol.map(calculateMinLST);  
  Export.table.toDrive({
    collection: imgCol, 
    description: "Drive-"+exportName,
    fileNamePrefix: exportName,
    folder: foldername,
    fileFormat: "CSV"
  });
  return imgCol;
}


function getDailyMaxDifferLST(start, finish, qTag, ndvi_threshold, MosaicFire, exportName,foldername) {
  var imgCol = getAllLSTImages(start, finish, qTag, ndvi_threshold, MosaicFire);       
  imgCol = imgCol.map(calculateMaxDifferLST);  
  Export.table.toDrive({
    collection: imgCol, 
    description: "Drive-"+exportName,
    fileNamePrefix: exportName,
    folder: foldername,
    fileFormat: "CSV"
  });
  return imgCol;
}


function getGlobalMaxLSTImage(start, finish, qTag, MosaicFire, exportName,foldername) {
  var imgCol = getDayLSTImages(start, finish, qTag, null, MosaicFire);
  var image = imgCol.reduce(ee.Reducer.max(2))
                    .select(["max", "max1"], ["lstmax", "lstmax_qc"]);
  image = image.addBands(ee.Image.pixelLonLat().select(["longitude","latitude"], ["lon", "lat"]));
  image = image.toDouble();
  Export.image.toDrive({
    image: image,
    description: "Drive-"+exportName,
    folder: foldername,
    fileNamePrefix: exportName,
    scale: 1000,
    region: roi,
    maxPixels: 1e13
  });
  return image;
}


function getGlobalMinLSTImage(start, finish, qTag, exportName,foldername) {
  var imgCol = getNightLSTImages(start, finish, qTag);
  var image = imgCol.reduce(ee.Reducer.min(2))
                    .select(["min", "min1"], ["lstmin", "lstmin_qc"]);
  image = image.addBands(ee.Image.pixelLonLat().select(["longitude","latitude"], ["lon", "lat"]));
  image = image.toDouble();
  Export.image.toDrive({
    image: image,
    description: "Drive-"+exportName,
    folder: foldername,
    fileNamePrefix: exportName,
    scale: 1000,
    region: roi,
    maxPixels: 1e13
  });
  return image;
}



function calculateMaxDifferWithinPeriod(image) {
  var qc_day = image.select("QC_Day");
  var qc_night = image.select("QC_Night");
  var lst_day = image.select("LST_Day_1km");
  var lst_night = image.select("LST_Night_1km");
  var differImage = lst_day.subtract(lst_night)
                          .abs()
                          .rename("differ");
  differImage = differImage.addBands(qc_day)
                           .addBands(qc_night);
  return differImage;
}



function getGlobalMaxDifLSTImage(start, finish, qTag, MosaicFire, ndvi_threshold, exportName,foldername) {
  // var imgColDay = getDayLSTImages(start, finish, qTag, null, MosaicFire);
  // var imgColNight = getNightLSTImages(start, finish, qTag);
  var imgCol = getAllLSTImages(start, finish, qTag, ndvi_threshold, MosaicFire);    
  imgCol = imgCol.map(calculateMaxDifferWithinPeriod); 
  // var image = imgCol.reduce(ee.Reducer.max(2))
  //                   .select(["differ", "differ1"], ["lstmax", "lstmax_qc"]);
  var image = imgCol.reduce(ee.Reducer.max(3))
                    .select(["max", "max1", "max2"], ["differ", "qc_day", "qc_night"]);
  image = image.addBands(ee.Image.pixelLonLat().select(["longitude","latitude"], ["lon", "lat"]));
  image = image.toDouble();
  Export.image.toDrive({
    image: image,
    description: "Drive-"+exportName,
    folder: foldername,
    fileNamePrefix: exportName,
    scale: 1000,
    region: roi,
    maxPixels: 1e13
  });
  return image;
}


function calculateAverageLST(image) {
  // var qc_day = image.select("QC_Day");
  // var qc_night = image.select("QC_Night");
  var lst_day = image.select("LST_Day_1km");
  var lst_night = image.select("LST_Night_1km");
  var averageImage = lst_day.add(lst_night)
                          .divide(2)
                          .rename("average");
                          
  var dict = averageImage.reduceRegion({
    reducer: ee.Reducer.mean()
               .setOutputs(["average"]),
    geometry: roi, 
    scale: 1000,
    maxPixels: 1e13
  });
  image = image.set("average", ee.Number(dict.get("average")));
  
  return image;
}

function getAverage(start, finish, qTag, ndvi_threshold, MosaicFire, exportName,foldername) {
  var imgCol = getAllLSTImages(start, finish, qTag, ndvi_threshold, MosaicFire);       
  imgCol = imgCol.map(calculateAverageLST);
  Export.table.toDrive({
    collection: imgCol, 
    description: "Drive-"+exportName,
    fileNamePrefix: exportName,
    folder: foldername,
    fileFormat: "CSV"
  });
  return imgCol;
}




/////////////////////////////////////////////////////////////

function main() {
  Map.addLayer(roi, {}, "roi", false);
  var startDateStr = "2002-7-1";
  var endDateStr = "2019-12-28";
  // var startDateStr = "2002-7-1";
  // var endDateStr = "2019-1-1";
  var start = ee.Date(startDateStr);
  var finish = ee.Date(endDateStr);
  //quality control
  var qTag = 1; //0 
 
  //NDVI and fire
  // var ndvi_threshold = 1;
  var ndvi_threshold = null;
  var MosaicFire = getMosaicFire(start, finish);
  //IGBP
  // var land_covers = [1, 16];
  // var MosaicIGBP = getMosaicIGBP(land_covers, start, finish);
  
  // export csv
  var exportName = "maxLST_"+startDateStr+"_"+endDateStr;
  var foldername = "MosaicMydMeanGlobal5";
  var maxLSTImgCol = getDailyMaxLST(start, finish, qTag, ndvi_threshold, MosaicFire, exportName, foldername);
  exportName = "minLST_"+startDateStr+"_"+endDateStr;
  var minLSTImgCol = getDailyMinLST(start, finish, qTag, exportName, foldername);
  exportName = "maxDifferLST_"+startDateStr+"_"+endDateStr;
  var maxDifferLSTImgCol = getDailyMaxDifferLST(start, finish, qTag, ndvi_threshold, MosaicFire, exportName, foldername);
  exportName = "meanLST_"+startDateStr+"_"+endDateStr;
  var meanLSTImage = getAverage(start, finish, qTag, ndvi_threshold, MosaicFire, exportName,foldername);
  // export image
  exportName = "maxLST_"+startDateStr+"_"+endDateStr;
  var maxLSTImage = getGlobalMaxLSTImage(start, finish, qTag, MosaicFire, exportName,foldername);
  exportName = "minLST_"+startDateStr+"_"+endDateStr;
  var minLSTImage = getGlobalMinLSTImage(start, finish, qTag, exportName,foldername);
  exportName = "maxdifLST_"+startDateStr+"_"+endDateStr;
  var maxDifLSTImage = getGlobalMaxDifLSTImage(start, finish, qTag, MosaicFire, ndvi_threshold, exportName,foldername);
  
  // // display layer to map
  // showGlobalLSTImage(maxLSTImage.select("lstmax"), "lstmax");
  // showGlobalLSTImage(minLSTImage.select("lstmin"), "lstmin");
  // showGlobalLSTImage(maxDifLSTImage.select("differ"), "differ");
  // showAndExportDailyMaxLSTPoint(maxLSTImgCol, false);
  // showAndExportDailyMinLSTPoint(minLSTImgCol, false);
  // showAndExportDailyMaxDifferLSTPoint(maxDifferLSTImgCol, false);
}

main();
