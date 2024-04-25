#%% Import libraries
import chunk
import sys
import xarray as xr
import numpy as np
import fiona
import geopandas as gp
import pandas as pd
from time import sleep
import matplotlib.pyplot as plt
import time
import dask
import shutil
from glob import glob
dask.config.set(**{'array.slicing.split_large_chunks': True})
import __utils

start_time = time.time()
chnk_sz = __utils.i_chnk_sz_space
chnk_time = __utils.i_chnk_sz

gage_id_attribute = __utils.gage_id_attribute_in_shapefile
#%% testing 
year = 2010
f_repo = "/scratch/dcl3nd/highres-radar-rainfall-processing/"
f_data = f_repo + "data/"
fpattern_in_ncs_mrms = f_data + "mrms_nc_preciprate_fullres_dailyfiles_constant_tstep/{}*.nc".format(year) # "${assar_dirs[out_fullres_dailyfiles]}${year}*.nc" 
fpattern_in_ncs_mrms_nonbiascorrected = f_data + "mrms_nc_preciprate_fullres_dailyfiles/{}*.nc".format(year) # "${assar_dirs[out_fullres_dailyfiles]}${year}*.nc" 
fpattern_in_stage_iv = "/scratch/dcl3nd/stormy/data/climate/StageIV_rainfall/{}/*.nc".format(year) # "${assar_dirs[stageiv_rainfall]}${year}/*.nc" | /scratch/dcl3nd/stormy/data/climate/StageIV_rainfall/2013/*.nc
shp_gages = f_data + "geospatial/rain_gages.shp" # ${assar_dirs[shp_gages]}
f_out_csv_mrms = f_data + "mrms_csv_preciprate_fullres_yearlyfiles_atgages/{}.csv".format(year) # "${assar_dirs[out_fullres_yearly_csvs_atgages]}${year}.csv"
f_out_nc_mrms = f_data + "mrms_nc_preciprate_fullres_yearlyfiles_atgages/{}.nc".format(year) # "${assar_dirs[out_fullres_yearly_atgages]}${year}.nc"
f_out_csv_stage_iv = f_data + "stage_iv_csv_preciprate_fullres_yearlyfiles_atgages/{}.csv".format(year) # "${assar_dirs[out_fullres_yearly_csvs_atgages_stageiv]}${year}.csv"
f_out_nc_stage_iv = f_data + "stage_iv_nc_preciprate_fullres_yearlyfiles_atgages/{}.nc".format(year) # "${assar_dirs[out_fullres_yearly_atgages_stageiv]}${year}.nc"

      
#%% load input parameters
fpattern_in_ncs_mrms = str(sys.argv[1])
fpattern_in_ncs_mrms_nonbiascorrected = str(sys.argv[2])
fpattern_in_stage_iv = str(sys.argv[3])
shp_gages = str(sys.argv[4])
f_out_csv_mrms = str(sys.argv[5])
f_out_nc_mrms = str(sys.argv[6])
f_out_csv_stage_iv = str(sys.argv[7])
f_out_nc_stage_iv = str(sys.argv[8])
#%% load data and define functions
gdf_gages = gp.read_file(shp_gages)

def find_closest_lat_lon(ar_lats, ar_lons, comp_lat, comp_lon):
    idx_closest_lat = np.argmin(abs(ar_lats - comp_lat))
    idx_closest_lon = np.argmin(abs(ar_lons - comp_lon))
    return idx_closest_lat, idx_closest_lon
#%% process mrms data
f_in_ncs_mrms = glob(fpattern_in_ncs_mrms)
f_in_ncs_mrms.sort()

f_in_ncs_mrms_uncrctd = glob(fpattern_in_ncs_mrms_nonbiascorrected)
f_in_ncs_mrms_uncrctd.sort()

ds_uncrctd = xr.open_mfdataset(f_in_ncs_mrms_uncrctd,  concat_dim = "time",
                chunks={'latitude':chnk_sz, 'longitude':chnk_sz},
                combine = "nested", engine = 'h5netcdf', coords='minimal')

if len(glob(f_in_ncs_mrms)) > 0:
    mrms = xr.open_mfdataset(f_in_ncs_mrms,  concat_dim = "time",
                chunks={'latitude':chnk_sz, 'longitude':chnk_sz},
                combine = "nested", engine = 'h5netcdf', coords='minimal')
    # print("Loaded data for {}. Time elapsed: {}".format(f_in_ncs, time.time()-start_time))
    # extract mrms data at gage locations
    # extract lat and long values and shift them from upper left reprsentation to
    # center of the rectangle
    if ds_uncrctd.attrs["gridcell_feature_represented_by_coordinate"] != "upper_left":
        sys.exit("User defined error. Script failed. The coordinates do not represent the upper left of each gridcell.")
    elif ds_uncrctd.attrs["longitude_units"] != "degrees_east":
        sys.exit("User defined error. Script failed. Longitude units are not degrees east.")
    elif ds_uncrctd.attrs["latitude_units"] != "degrees_north":
        sys.exit("User defined error. Script failed. Latitude units are not degrees north.")
    else:
        pass
    # shift gridcells from upper left to center
    ## MRMS
    mrms_lat_upper_left = mrms["latitude"].values
    mrms_lat_centered = mrms_lat_upper_left - ds_uncrctd.attrs['grid_spacing']/2 # shift coordinates south 1/2 a gridcell
    mrms_long_upper_left = mrms["longitude"].values
    mrms_long_centered = mrms_long_upper_left + ds_uncrctd.attrs['grid_spacing']/2 # shift coordaintes east 1/2 a gridcell
    # add lat and long attribute to gage geodataframe to match with mrms cells
    for i in range(0, len(gdf_gages)):
        xy = list(gdf_gages.to_crs("EPSG:4326").iloc[i, :].geometry.coords)
        long = gdf_gages.loc[i, "long"] = xy[0][0] + 360 # convert from negative degrees east to positive degrees east
        lat = gdf_gages.loc[i, "lat"] = xy[0][1]
        # finding the index is based on the centered lat and long
        # mrms_idx_closest_lat = np.argmin(abs(mrms_lat_centered - lat))
        # mrms_idx_closest_lon = np.argmin(abs(mrms_long_centered - long))
        mrms_idx_closest_lat, mrms_idx_closest_lon = find_closest_lat_lon(mrms_lat_centered, mrms_long_centered, lat, long)
        # assign gdg_gages attributes
        gdf_gages.loc[i, "mrms_lat_idx"] = mrms_idx_closest_lat
        gdf_gages.loc[i, "mrms_long_idx"] = mrms_idx_closest_lon
        # this should be the original lat and lon in the netcdf (representing upper left)
        gdf_gages.loc[i, "mrms_lat"] = mrms_lat_upper_left[mrms_idx_closest_lat]
        gdf_gages.loc[i, "mrms_long"] = mrms_long_upper_left[mrms_idx_closest_lon]

    mrms_lat_idx_min = min(gdf_gages['mrms_lat_idx'])
    mrms_lat_idx_max = max(gdf_gages['mrms_lat_idx'])
    mrms_long_idx_min = min(gdf_gages['mrms_long_idx'])
    mrms_long_idx_max = max(gdf_gages['mrms_long_idx'])

    # subset mrms data with a 5 gridcell buffer
    mrms_lat_idx = np.arange(mrms_lat_idx_min-5, mrms_lat_idx_max+6, dtype=int) 
    mrms_long_idx = np.arange(mrms_long_idx_min-5, mrms_long_idx_max+6, dtype=int)


    idx_mrms = dict(latitude = mrms_lat_idx, longitude = mrms_long_idx)
    event_data_mrms = mrms[idx_mrms]
    event_data_mrms = event_data_mrms.chunk(chunks={'latitude':chnk_sz, 'longitude':chnk_sz})

    # print("extracting gage_id variable...")
    gage_ids_mrms = gdf_gages.loc[:, ['mrms_lat', 'mrms_long', gage_id_attribute]]

    # print("Loading event_data into memory...")
    bm_time = time.time()
    event_data_mrms_loaded = event_data_mrms.load()
    # print("Time to load subset of data into memory (s): {}".format(time.time() - bm_time))
    # Process dataframe and export to csv
    # print("Creating csv...")
    bm_time = time.time()

    df_mrms = event_data_mrms_loaded.rainrate.to_dataframe()

    df_mrms_rst_ind = df_mrms.reset_index()


    df_out_mrms = pd.merge(df_mrms_rst_ind, gage_ids_mrms, how = 'left', left_on = ['latitude', 'longitude'], right_on = ['mrms_lat', 'mrms_long'])

    # drop columns that don't coincide with a rain gage
    df_out_mrms.dropna(subset=[gage_id_attribute], inplace=True)
    df_out_mrms = df_out_mrms.reset_index()
    df_out_mrms.rename(columns={"rainrate":"precip_mm_per_hour", gage_id_attribute:"overlapping_gage_id"}, inplace = True)
    df_out_mrms = df_out_mrms.loc[:, ['time', 'precip_mm_per_hour', 'mrms_lat', 'mrms_long', 'overlapping_gage_id']]

    # export
    df_out_mrms.to_csv(f_out_csv_mrms)

    # print("Time to export csv (s): {}".format(time.time() - bm_time))
    # export netcdf
    # print("Creating netcdf...")
    bm_time = time.time()
    event_data_mrms_loaded.to_netcdf(f_out_nc_mrms)

else:
    print("no mrms data for this year")
#%% process stage IV
f_in_stage_iv = glob(fpattern_in_stage_iv)
f_in_stage_iv.sort()
if len(glob(f_in_stage_iv)) > 0:
    try:
        stage_iv = xr.open_mfdataset(f_in_stage_iv,  concat_dim = "time",
                    chunks={'outlat':chnk_sz, 'outlon':chnk_sz},
                    combine = "nested", coords='minimal', engine = 'h5netcdf')
    except:
        print("loading stage IV data failed for some reason.....")
        print("files attempted to be loaded: ")
        print(glob(f_in_stage_iv))

    # format stage iv 
    new_lon = stage_iv["longitude"].values[0,:]+360 # convert from degrees west to degrees east
    new_lat = stage_iv["latitude"].values[0,:]
    stage_iv['outlon'] =  new_lon
    stage_iv['outlat'] =  new_lat
    stage_iv = stage_iv.drop_vars(["latitude", "longitude"])
    stage_iv = stage_iv.rename_dims(dims_dict=dict(outlat = "latitude", outlon="longitude"))
    stage_iv = stage_iv.rename(dict(outlat = "latitude", outlon = "longitude"))

    # define grid spacing attribute
    lat_diff = np.unique(np.diff(stage_iv.latitude.values))
    lon_diff = np.unique(np.diff(stage_iv.longitude.values))

    # if abs(lat_diff) != abs(lon_diff):
    #     print("stage_iv.latitude.values")
    #     print(stage_iv.latitude.values)
    #     print("##################################################")
    #     print("stage_iv.latitude.values")
    #     print(stage_iv.latitude.values)
    #     print("##################################################")
    #     print("lat_diff")
    #     print(lat_diff)
    #     print("##################################################")
    #     print("lon_diff")
    #     print(lon_diff)
    #     print("warning: Stage IV latitude grid spacing is irregular does not align with longitude grid spacing")

    stage_iv.latitude.attrs['units'] = 'degrees_north'
    stage_iv.longitude.attrs['units'] = 'degrees_east'
    stage_iv.attrs['grid_spacing'] = float(abs(lat_diff)[0])

    ## stage IV
    stage_iv_lat_upper_left = stage_iv["latitude"].values
    stage_iv_lat_centered = stage_iv_lat_upper_left - stage_iv.attrs['grid_spacing']/2 # shift coordinates south 1/2 a gridcell
    stage_iv_long_upper_left = stage_iv["longitude"].values
    stage_iv_long_centered = stage_iv_long_upper_left + stage_iv.attrs['grid_spacing']/2 # shift coordaintes east 1/2 a gridcell

    # add lat and long attribute to gage geodataframe to match with mrms cells
    for i in range(0, len(gdf_gages)):
        xy = list(gdf_gages.to_crs("EPSG:4326").iloc[i, :].geometry.coords)
        long = gdf_gages.loc[i, "long"] = xy[0][0] + 360 # convert from negative degrees east to positive degrees east
        lat = gdf_gages.loc[i, "lat"] = xy[0][1]
        stage_iv_idx_closest_lat, stage_iv_idx_closest_lon = find_closest_lat_lon(stage_iv_lat_centered, stage_iv_long_centered, lat, long)
        gdf_gages.loc[i, "stage_iv_lat_idx"] = stage_iv_idx_closest_lat
        gdf_gages.loc[i, "stage_iv_long_idx"] = stage_iv_idx_closest_lon
        gdf_gages.loc[i, "stage_iv_lat"] = stage_iv_lat_upper_left[stage_iv_idx_closest_lat]
        gdf_gages.loc[i, "stage_iv_long"] = stage_iv_long_upper_left[stage_iv_idx_closest_lon]

    stage_iv_lat_idx_min = min(gdf_gages['stage_iv_lat_idx'])
    stage_iv_lat_idx_max = max(gdf_gages['stage_iv_lat_idx'])
    stage_iv_long_idx_min = min(gdf_gages['stage_iv_long_idx'])
    stage_iv_long_idx_max = max(gdf_gages['stage_iv_long_idx'])

    # stage iv
    stage_iv_lat_idx = np.arange(stage_iv_lat_idx_min-5, stage_iv_lat_idx_max+6, dtype=int) 
    stage_iv_long_idx = np.arange(stage_iv_long_idx_min-5, stage_iv_long_idx_max+6, dtype=int)

    idx_stage_iv = dict(latitude = stage_iv_lat_idx, longitude = stage_iv_long_idx)
    event_data_stage_iv = stage_iv[idx_stage_iv]
    event_data_stage_iv = event_data_stage_iv.chunk(chunks={'latitude':chnk_sz, 'longitude':chnk_sz})

    gage_ids_stage_iv = gdf_gages.loc[:, ['stage_iv_lat', 'stage_iv_long', gage_id_attribute]]

    event_data_stage_iv_loaded = event_data_stage_iv.load()

    df_stage_iv = event_data_stage_iv_loaded.rainrate.to_dataframe(dim_order = ["latitude", "longitude", "time"])
    # df_stage_iv = pd.DataFrame(dict(latitude = event_data_stage_iv_loaded.latitude.values,
    #                                 longitude = event_data_stage_iv_loaded.longitude.values,
    #                                 rainrate = event_data_stage_iv_loaded.rainrage.values))

    df_stage_iv_rst_ind = df_stage_iv.reset_index()
    # replace lat and lon indices with actual values
    lats = event_data_stage_iv_loaded.latitude.values
    lons = event_data_stage_iv_loaded.longitude.values

    lat_col = df_stage_iv_rst_ind.latitude.replace(df_stage_iv_rst_ind.latitude.unique(), lats)
    lon_col = df_stage_iv_rst_ind.longitude.replace(df_stage_iv_rst_ind.longitude.unique(), lons)

    df_stage_iv_rst_ind["latitude"] = lat_col
    df_stage_iv_rst_ind["longitude"] = lon_col

    df_out_stage_iv = pd.merge(df_stage_iv_rst_ind, gage_ids_stage_iv, how = 'left', left_on = ['latitude', 'longitude'], right_on = ['stage_iv_lat', 'stage_iv_long'])

    df_out_stage_iv.dropna(subset=[gage_id_attribute], inplace=True)
    df_out_stage_iv = df_out_stage_iv.reset_index()
    df_out_stage_iv.rename(columns={"rainrate":"precip_mm_per_hour", gage_id_attribute:"overlapping_gage_id"}, inplace = True)
    df_out_stage_iv = df_out_stage_iv.loc[:, ['time', 'precip_mm_per_hour', 'stage_iv_lat', 'stage_iv_long', 'overlapping_gage_id']]

    # export
    df_out_stage_iv.to_csv(f_out_csv_stage_iv)

    event_data_stage_iv_loaded.to_netcdf(f_out_nc_stage_iv)
else:
    print("no stage iv data for this year")
#%% performance
# print("Time to export netcdf (s): {}".format(time.time() - bm_time))
print("Finished creating {}, {}, {}, and {}. Time elapsed: {}".format(f_out_nc_mrms, f_out_nc_stage_iv, f_out_csv_mrms, f_out_csv_stage_iv, time.time()-start_time))