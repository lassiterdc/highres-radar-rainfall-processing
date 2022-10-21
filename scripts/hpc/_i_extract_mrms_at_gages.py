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
dask.config.set(**{'array.slicing.split_large_chunks': True})
import __utils

start_time = time.time()
chnk_sz = __utils.i_chnk_sz_space
chnk_time = __utils.i_chnk_sz

gage_id_attribute = __utils.gage_id_attribute_in_shapefile
#%% load input parameters
f_in_ncs = str(sys.argv[1])
shp_gages = str(sys.argv[2])
f_out_csv = str(sys.argv[3])
f_out_nc = str(sys.argv[4])
#%% load data
gdf_gages = gp.read_file(shp_gages)

mrms = xr.open_mfdataset(f_in_ncs,  concat_dim = "time",
            chunks={'latitude':chnk_sz, 'longitude':chnk_sz},
            combine = "nested", engine = 'h5netcdf', coords='minimal')

# print("Loaded data for {}. Time elapsed: {}".format(f_in_ncs, time.time()-start_time))
#%% extract mrms data at gage locations
# extract lat and long values and shift them from upper left reprsentation to
# center of the rectangle
if mrms.attrs["gridcell_feature_represented_by_coordinate"] != "upper_left":
    sys.exit("User defined error. Script failed. The coordinates do not represent the upper left of each gridcell.")

if mrms.attrs["longitude_units"] != "degrees_east":
    sys.exit("User defined error. Script failed. Longitude units are not degrees east.")

if mrms.attrs["latitude_units"] != "degrees_north":
    sys.exit("User defined error. Script failed. Latitude units are not degrees north.")

mrms_lat_upper_left = mrms["latitude"].values
mrms_lat_centered = mrms_lat_upper_left - mrms.attrs['grid_spacing']/2 # shift coordinates south 1/2 a gridcell
mrms_long_upper_left = mrms["longitude"].values
mrms_long_centered = mrms_long_upper_left + mrms.attrs['grid_spacing']/2 # shift coordaintes east 1/2 a gridcell

# add lat and long attribute to gage geodataframe to match with mrms cells
for i in range(0, len(gdf_gages)):
    xy = list(gdf_gages.to_crs("EPSG:4326").iloc[i, :].geometry.coords)
    long = gdf_gages.loc[i, "long"] = xy[0][0] + 360 # convert from negative degrees east to positive degrees east
    lat = gdf_gages.loc[i, "lat"] = xy[0][1]
    # finding the index is based on the centered lat and long
    mrms_idx_closest_lat = np.argmin(abs(mrms_lat_centered - lat))
    mrms_idx_closest_lon = np.argmin(abs(mrms_long_centered - long))
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

idx = dict(latitude = mrms_lat_idx, longitude = mrms_long_idx)

event_data = mrms[idx]
event_data = event_data.chunk(chunks={'latitude':chnk_sz, 'longitude':chnk_sz})

# print("extracting gage_id variable...")
gage_ids = gdf_gages.loc[:, ['mrms_lat', 'mrms_long', gage_id_attribute]]

# print("Loading event_data into memory...")
bm_time = time.time()
event_data_loaded = event_data.load()
print("Time to load subset of data into memory (s): {}".format(time.time() - bm_time))
#%% Process dataframe and export to csv
print("Creating csv...")
bm_time = time.time()

df = event_data_loaded.rainrate.to_dataframe()

df_rst_ind = df.reset_index()

df_out = pd.merge(df_rst_ind, gage_ids, how = 'left', left_on = ['latitude', 'longitude'], right_on = ['mrms_lat', 'mrms_long'])

# drop columns that don't coincide with a rain gage
df_out.dropna(subset=[gage_id_attribute], inplace=True)
df_out = df_out.reset_index()
df_out.rename(columns={"rainrate":"precip_mm_per_hour", gage_id_attribute:"overlapping_gage_id"}, inplace = True)
df_out = df_out.loc[:, ['time', 'precip_mm_per_hour', 'mrms_lat', 'mrms_long', 'overlapping_gage_id']]
df_out.to_csv(f_out_csv)

# print("Time to export csv (s): {}".format(time.time() - bm_time))
#%% export netcdf
# print("Creating netcdf...")
bm_time = time.time()
event_data_loaded.to_netcdf(f_out_nc)
# print("Time to export netcdf (s): {}".format(time.time() - bm_time))
print("Finished creating csv and nc for. Time elapsed: {}".format(time.time()-start_time))