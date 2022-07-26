#!/usr/bin/env python

### Import Modules ###

import sys
import os
import geopandas as gpd
import numpy as np
from osgeo import gdal
from osgeo import ogr
from osgeo import gdalconst
import fiona
import rasterio
import rasterio.mask
import rasterio
import rasterio.mask

# Pass run detail

path_local = sys.argv[1]
path_run = sys.argv[2]
cm = int(sys.argv[3])


#%% 

def mkdir(dir):
    
    if not os.path.exists(dir):
        os.mkdir(dir)

path_data = path_run + 'data/'

path_RICorDE_processed = path_data + 'RICorDE_processed/'
path_RICorDE_processed_layer = path_RICorDE_processed + 'RICorDE_layer/'
path_RICorDE_processed_layer_tif_processed = path_RICorDE_processed_layer + 'tif/'
path_RICorDE_processed_layer_shp_processed = path_RICorDE_processed_layer + 'shp/'

crs = {'init' :'epsg:4326'} 

#%%
  
cm = cm/100
print(cm)

RICorDE_shape = path_RICorDE_processed_layer_shp_processed + f"idai_r1_0331_depths_resampled_9as_reproj_{cm}.shp"

mask = None

from rasterio.features import shapes

print(path_RICorDE_processed + "idai_r1_0331_depths_resampled_9as_reproj.tif")

with rasterio.Env():

    with rasterio.open(path_RICorDE_processed + "idai_r1_0331_depths_resampled_9as_reproj.tif") as src:
        image = src.read()
        image[np.isnan(image)] = 0
        image[image<cm] = 0
        image[image>=cm] = 1


    results = (
    {'properties': {'raster_val': v}, 'geometry': s}
    for i, (s, v)
    in enumerate(
    shapes(image, mask=mask, transform=src.transform)))

    with fiona.open(
    RICorDE_shape, 'w',
    driver='Shapefile',
    crs=src.crs,
    schema={'properties': [('raster_val', 'int')],
    'geometry': 'Polygon'}) as dst:
        dst.writerecords(results)

satellite_shape_dissolved = gpd.read_file(RICorDE_shape)
satellite_shape_dissolved = satellite_shape_dissolved[satellite_shape_dissolved['raster_val'] == 1]
satellite_shape_dissolved.to_file(RICorDE_shape)

rst_fn = path_RICorDE_processed + '/cfwindzos065_no_dif_cropped_2.tif'
shp_fn = RICorDE_shape

#This raster is the model for our output (CRS, extent)
ndsm = rst_fn
#This shapefile contains the features we want to burn
shp = shp_fn

data = gdal.Open(ndsm, gdalconst.GA_ReadOnly)

geo_transform = data.GetGeoTransform()
x_min = geo_transform[0]
y_max = geo_transform[3]
x_max = x_min + geo_transform[1] * data.RasterXSize
y_min = y_max + geo_transform[5] * data.RasterYSize
x_res = data.RasterXSize
y_res = data.RasterYSize


vec = ogr.Open(shp)
lyr = vec.GetLayer(0)
pixel_width = geo_transform[1]

output = path_RICorDE_processed_layer_tif_processed + f"idai_r1_0331_depths_resampled_9as_reproj_{cm}.tif"
target_ds = gdal.GetDriverByName('GTiff').Create(output, x_res, y_res, 1, gdal.GDT_UInt32)
target_ds.SetGeoTransform((x_min, pixel_width, 0, y_min, 0, pixel_width))
band = target_ds.GetRasterBand(1)
NoData_value = 0
band.SetNoDataValue(NoData_value)
band.FlushCache()


driver = ogr.GetDriverByName('ESRI Shapefile')
for feat in lyr:
    burn_value = int(feat.GetField("raster_val"))
    datasource = driver.CreateDataSource(f"{path_data}RICorDE_processed/RICorDE_layer/shp_temp/{cm}_temp.shp")
    layer = datasource.CreateLayer('temp', lyr.GetSpatialRef(), geom_type=ogr.wkbPolygon)
    layer.CreateFeature(feat)
    gdal.RasterizeLayer(target_ds, [1], layer, burn_values=[burn_value], options=["ALL_TOUCHED=FALSE"])
    datasource.Destroy()

#This makes raster to write to disk
target_ds = None    

#%%

print('job succesfully completed!')

