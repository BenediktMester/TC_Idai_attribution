#!/usr/bin/env python

### Import Modules ###

import sys
import os
import geopandas as gpd
import pandas as pd
from rasterstats import zonal_stats
import numpy as np
import rasterio
import fiona

# Pass run detail

path_local = sys.argv[1]
path_run = sys.argv[2]
tide = sys.argv[3]
flood_level_threshold = sys.argv[4]

#%%
 
def mkdir(dir):
    
    if not os.path.exists(dir):
        os.mkdir(dir)     

path_data = path_run + 'data/'

path_data_processed = path_run + 'data_processed/'
path_data_processed_creator = mkdir(path_data_processed)
path_GeoClaw_processed = path_data_processed + 'GeoClaw_processed/'
path_GeoClaw_processed_creator = mkdir(path_GeoClaw_processed)
path_GeoClaw_factual_cropped = path_data_processed + 'GeoClaw_processed/GeoClaw_factual_cropped/'
path_GeoClaw_factual_cropped_creator = mkdir(path_GeoClaw_factual_cropped)
path_flood_product_processed = path_data_processed + 'flood_product_processed/'
path_flood_product_processed_creator = mkdir(path_flood_product_processed)
path_flood_product_merged = path_flood_product_processed + 'flood_product_merged/'
path_flood_product_merged_creator = mkdir(path_flood_product_merged)
path_flood_product_merged_bi = path_flood_product_processed + 'flood_product_merged_bi/'
path_flood_product_merged_creator = mkdir(path_flood_product_merged_bi)
path_population_processed = path_data_processed + 'population_processed/'
path_population_processed_creator = mkdir(path_population_processed)
path_results = path_run + 'results/'
path_results_creator = mkdir(path_results)

crs = {'init' :'epsg:4326'} 

#%%

profile = rasterio.open(path_data + 'GeoClaw/factual/2019063S18038_cf-zos_aviso-fes_mean.tif').profile

for filename in sorted(os.listdir(path_data + 'GeoClaw/factual/')):

    if tide in filename:

        print(filename)
        
        # Clip Factual
        
        filename_f = filename.split('.tif')[0]
                
        mask = None  

        GeoClaw_shape = path_data + '/GeoClaw_shp_template/GeoClaw_clipped.shp'

        with fiona.open(GeoClaw_shape, "r") as shapefile:
            shapes = [feature["geometry"] for feature in shapefile]

        with rasterio.open(f"{path_data}GeoClaw/factual/{filename_f}.tif") as src:
            out_image, out_transform = rasterio.mask.mask(src, shapes, crop=True)
            out_meta = src.meta

        out_meta.update({"driver": "GTiff",
                          "height": out_image.shape[1],
                          "width": out_image.shape[2],
                          "transform": out_transform})

        out_image[out_image==0] = 0 
        coast_f_dif_cropped = out_image.copy()

        with rasterio.open(f"{path_GeoClaw_factual_cropped}{filename_f}_cropped.tif", "w", **out_meta) as dest:
                dest.write(out_image) 
                
for filename in sorted(os.listdir(path_data + 'GeoClaw/counterfactual/')):

    if tide in filename:

        coast_f_count_temp = rasterio.open(path_data + 'GeoClaw/counterfactual/' + filename).read()
        coast_f_count_temp = np.round(coast_f_count_temp,1)

        filename_cut = filename.split('8038_')[1]
        print(filename_cut)

        filename_cut_cf = filename_cut.split('-')[0]
        filename_cut_tide = filename_cut.split('_')[2]
        filename_cut_tide = filename_cut_tide.split('.')[0]
   
        coast_f_fact = rasterio.open(f"{path_data}GeoClaw/factual/{filename_f}.tif").read()
        coast_f_fact = np.round(coast_f_fact,1)
        coast_f_dif = coast_f_fact - coast_f_count_temp 

        with rasterio.open(path_GeoClaw_processed + filename_cut_cf + '_' + filename_cut_tide + '_dif.tif', "w", **profile) as dest:
            dest.write(coast_f_dif) 

        # Clip GeoClaw

        mask = None  

        path_raster_layer_resampled = path_GeoClaw_processed + filename_cut_cf + '_' + filename_cut_tide + '_dif.tif'
        GeoClaw_shape = path_data + '/GeoClaw_shp_template/GeoClaw_clipped.shp'

        with fiona.open(GeoClaw_shape, "r") as shapefile:
            shapes = [feature["geometry"] for feature in shapefile]

        with rasterio.open(path_raster_layer_resampled) as src:
            out_image, out_transform = rasterio.mask.mask(src, shapes, crop=True)
            out_meta = src.meta

        out_meta.update({"driver": "GTiff",
                          "height": out_image.shape[1],
                          "width": out_image.shape[2],
                          "transform": out_transform})

        out_image[out_image==0] = 0 
        coast_f_dif_cropped = out_image.copy()

        path_raster_layer_resampled_cropped = path_GeoClaw_processed + filename_cut_cf + '_' + filename_cut_tide + '_dif_cropped.tif'

        with rasterio.open(path_raster_layer_resampled_cropped, "w", **out_meta) as dest:
                dest.write(out_image) 

        # Subtract satellite imagery

        satellite = rasterio.open(f"{path_data}RICorDE/idai_r1_0331_depths_0_resampled_9as_reproj_final.tif").read()

        satellite = np.fliplr(np.flipud(satellite))
        satellite[np.isnan(satellite)] = 0
        
        coast_f_fact = rasterio.open(f"{path_GeoClaw_factual_cropped}{filename_f}_cropped.tif").read()
        coast_f_fact = np.round(coast_f_fact,1)

        flood_f = np.maximum(satellite,coast_f_fact)  

        satellite_cf = flood_f - coast_f_dif_cropped

        
        satellite_cf[satellite_cf < 0] = 0
        satellite_profile = rasterio.open(f"{path_data}RICorDE/idai_r1_0331_depths_0_resampled_9as_reproj_final.tif").profile

        with rasterio.open(path_flood_product_merged + 'flood_product_merged_' + filename_cut_cf + '_' + filename_cut_tide + '.tif', "w", **out_meta) as dest:
            dest.write(satellite_cf) 

        for flood_level_threshold_temp in ([10,50,100]):

            print(flood_level_threshold_temp)

            satellite_cf_bi = satellite_cf.copy()

            satellite_cf_bi[satellite_cf_bi < flood_level_threshold_temp/100] = 0
            satellite_cf_bi[satellite_cf_bi > 0] = 1

            with rasterio.open(path_flood_product_merged_bi + 'flood_product_merged_' + filename_cut_cf + '_' + filename_cut_tide + f"_{flood_level_threshold_temp}_bi.tif", "w", **out_meta) as dest:
                dest.write(satellite_cf_bi) 
                
            satellite_bi_cf_shp = path_flood_product_merged_bi + 'flood_product_merged_' + filename_cut_cf + '_' + filename_cut_tide + f"_{flood_level_threshold_temp}_bi.shp"   
            mask = None
            from rasterio.features import shapes

            with rasterio.Env():

                with rasterio.open(path_flood_product_merged_bi + 'flood_product_merged_' + filename_cut_cf + '_' + filename_cut_tide + f"_{flood_level_threshold_temp}_bi.tif") as src:
                    image = src.read()

                results = (
                {'properties': {'raster_val': v}, 'geometry': s}
                for i, (s, v)
                in enumerate(
                shapes(image, mask=mask, transform=src.transform)))

                with fiona.open(
                satellite_bi_cf_shp, 'w',
                driver='Shapefile',
                crs=src.crs,
                schema={'properties': [('raster_val', 'int')],
                'geometry': 'Polygon'}) as dst:
                    dst.writerecords(results)

            satellite_shape_dissolved = gpd.read_file(satellite_bi_cf_shp)
            satellite_shape_dissolved = satellite_shape_dissolved[satellite_shape_dissolved['raster_val'] == 1]
            satellite_shape_dissolved.to_file(satellite_bi_cf_shp)

#%%

# Affected people

# Prepare population data

countries = gpd.read_file(path_local + 'data/states/gadm36_levels_shp/gadm36_0.shp')
MOZ = countries[countries.GID_0 == 'MOZ']

MOZ_2015 = 27042001
MOZ_2019 = 30366043

MOZ_pop_growth_factor = MOZ_2019 / MOZ_2015
print('MOZ_pop_growth_factor:',MOZ_pop_growth_factor)

GeoClaw_shape_geometry = gpd.read_file(path_data + '/GeoClaw_shp_template/GeoClaw_clipped.shp')['geometry']

with rasterio.open(path_local + "data/population/GHSL/GHS_POP_E2015_GLOBE_R2019A_4326_9ss_V1_0.tif") as src:
    out_image, out_transform = rasterio.mask.mask(src, GeoClaw_shape_geometry, crop=True)
    out_meta = src.meta

out_meta.update({"driver": "GTiff",
  "height": out_image.shape[1],
  "width": out_image.shape[2],
  "transform": out_transform})

out_image = out_image * MOZ_pop_growth_factor
out_image[out_image < 0] = 0

path_pop_tif_cropped = path_population_processed + 'MOZ_GHSL_9as_clipped.tif'

with rasterio.open(path_pop_tif_cropped, "w", **out_meta) as dest:
    dest.write(out_image) 

#%%

path_pop_tif_cropped = path_population_processed + 'MOZ_GHSL_9as_clipped.tif'
  
scenario_list = []
tide_list = []
flood_level_threshold_list = []
affected_list = []
 
print('to match:', tide, flood_level_threshold)
for filename in sorted(os.listdir(path_flood_product_merged_bi)): 
    
    print(filename)
    
    print(filename.split('_')[4])
    
    print(filename.split('_')[5])#.split('.')[0])
    
    print(filename.split('.')[1])
    
    if filename.split('_')[4] == tide and filename.split('_')[5] == flood_level_threshold and filename.split('.')[1] == 'shp':
            
        filename_cut_cf = filename.split('_')[3]
        
        print('it is a match!')

        print(filename_cut_cf)

        stat = 0


        flood = gpd.read_file(path_flood_product_merged_bi + 'flood_product_merged_' + filename_cut_cf + f"_{tide}_{flood_level_threshold}_bi.shp")
        try:
            stat = zonal_stats(flood,
                path_pop_tif_cropped, 
                all_touched = False,
                stats = ['sum'])                           

            stat = int(pd.DataFrame(stat).sum().iloc[0])
        except:
            stat = 0
        print(stat)
        
        scenario_list.append(filename_cut_cf)
        tide_list.append(tide)
        flood_level_threshold_list.append(flood_level_threshold.split('.')[0])
        affected_list.append(stat)

        ##
        results = pd.DataFrame()
        results.loc[:,'scenario'] = scenario_list
        results.loc[:,'tide'] = tide_list
        results.loc[:,'threshold'] = flood_level_threshold_list
        results.loc[:,'affected'] = affected_list
        results.to_csv(f"{path_results}results_{tide}_{flood_level_threshold}_temp.csv", index=False)
        
results = pd.DataFrame()
results.loc[:,'scenario'] = scenario_list
results.loc[:,'tide'] = tide_list
results.loc[:,'threshold'] = flood_level_threshold_list
results.loc[:,'affected'] = affected_list
results.to_csv(f"{path_results}results_{tide}_{flood_level_threshold}.csv", index=False)

#%%

# affected by high wind speeds

results = pd.read_csv(f"{path_results}results_{tide}_{flood_level_threshold}.csv")
results['affected_wind_64'] = np.nan
results['affected_wind_64_cf'] = np.nan
results['affected_wind_96'] = np.nan
results['affected_wind_96_cf'] = np.nan

path_TC_Idai_windspeed_processed_no_flood = path_data_processed + 'TC_Idai_windspeed_processed_no_flood/'
path_TC_Idai_windspeed_processed_no_flood_creator = mkdir(path_TC_Idai_windspeed_processed_no_flood)

Idai_TC_shp_dissolved_64 = gpd.read_file(path_data + 'TC_Idai_windspeed_processed/Idai_TC_final_64.shp')
Idai_TC_shp_dissolved_64_cf = gpd.read_file(path_data + 'TC_Idai_windspeed_processed/Idai_TC_final_64_cf.shp')
Idai_TC_shp_dissolved_96 = gpd.read_file(path_data + 'TC_Idai_windspeed_processed/Idai_TC_final_96.shp')

for i in range(len(results)):
    
    scenario_temp = results.loc[i,'scenario']
    print(results.loc[i,'scenario'])
    print(tide)
    print(flood_level_threshold)
    
    for filename in sorted(os.listdir(path_flood_product_merged_bi)):

        if filename == f"flood_product_merged_{scenario_temp}_{tide}_{flood_level_threshold}_bi.shp": 
            
            print('yes')
            print(filename)
            print(f"flood_product_merged_{scenario_temp}_{tide}_{flood_level_threshold}_bi.shp")    
            print(path_flood_product_merged_bi + filename)
            flood_temp = gpd.read_file(path_flood_product_merged_bi + filename)
            print("len:",len(flood_temp))
            print("len:",len(Idai_TC_shp_dissolved_64))

            TC_windspeed_no_flood_64 = gpd.overlay(Idai_TC_shp_dissolved_64,flood_temp,how='difference')
            TC_windspeed_no_flood_64.to_file(path_TC_Idai_windspeed_processed_no_flood + f"TC_2019_max_cropped_binary_final_64_no_flood_{scenario_temp}_{tide}_{flood_level_threshold}_bi.shp")

            TC_windspeed_no_flood_64_cf = gpd.overlay(Idai_TC_shp_dissolved_64_cf,flood_temp,how='difference')
            TC_windspeed_no_flood_64_cf.to_file(path_TC_Idai_windspeed_processed_no_flood + f"TC_2019_max_cropped_binary_final_64_cf_no_flood_{scenario_temp}_{tide}_{flood_level_threshold}_bi.shp")

            TC_windspeed_no_flood_96 = gpd.overlay(Idai_TC_shp_dissolved_96,flood_temp,how='difference')
            TC_windspeed_no_flood_96.to_file(path_TC_Idai_windspeed_processed_no_flood + f"TC_2019_max_cropped_binary_final_96_no_flood_{scenario_temp}_{tide}_{flood_level_threshold}_bi.shp")

            path_pop_tif_cropped = path_data_processed + "population_processed/MOZ_GHSL_9as_clipped.tif"

            # affected people by TC windspeeds

            try:
                stat = zonal_stats(TC_windspeed_no_flood_64,
                    path_pop_tif_cropped, 
                    all_touched = False,
                    stats = ['sum'])                           

                stat = int(pd.DataFrame(stat).sum().iloc[0])
            except:
                stat = 0
            print(stat)
            
            results.loc[i,'affected_wind_64'] = stat

            # affected people by TC windspeeds - cf
            try:
                stat = zonal_stats(TC_windspeed_no_flood_64_cf,
                    path_pop_tif_cropped, 
                    all_touched = False,
                    stats = ['sum'])                           

                stat = int(pd.DataFrame(stat).sum().iloc[0])
            except:
                stat = 0
            print(stat)
            
            results.loc[i,'affected_wind_64_cf'] = stat
            
            # affected people by TC windspeeds 96 kn

            try:
                stat = zonal_stats(TC_windspeed_no_flood_96,
                    path_pop_tif_cropped, 
                    all_touched = False,
                    stats = ['sum'])                           

                stat = int(pd.DataFrame(stat).sum().iloc[0])
            except:
                stat = 0
            print(stat)
            
            results.loc[i,'affected_wind_96'] = stat
            results.loc[i,'affected_wind_96_cf'] = 0

            results.to_csv(f"{path_results}results_{tide}_{flood_level_threshold}_wind_temp.csv")

#%%

results.to_csv(f"{path_results}results_{tide}_{flood_level_threshold}_wind_final.csv", index=False)

#%%

print('job succesfully completed!')