"""
Goal is to create hourly and daily versions of the dataset so it's easy to compare to other data sets
"""

#%% libraries and directories
import time
start_time = time.time()
from doctest import debug_script
from importlib.util import decode_source
import xarray as xr
from glob import glob
import pandas as pd
import numpy as np
import dask
dask.config.set(**{'array.slicing.split_large_chunks': False})
import shutil
import sys

#%% parameters
chnk_sz = "5000MB"
#%% testing on local pc
# in_date = 20210503
# f_in_nc = "data/mrms_for_rainyday_subset/" + "{}.nc".format(in_date)
# fl_out_zar = "out_zarr/" + "_q_{}_hourly.zarr".format(in_date)
# f_out_nc_hrly = "out_netcdfs/_q_{}_hourly.nc".format(in_date)
# f_out_nc_daily = "out_netcdfs/_q_{}_daily.nc".format(in_date)
# chnk_sz = "5000MB"

#%% testing on rivanna interactive job from command line
# module purge
# module load anaconda
# source activate mrms_processing
# python _d4a_resampling_to_hourly_and_daily_timesteps.py 20210504 "/project/quinnlab/dcl3nd/norfolk/data/mrms_for_rainyday/" "/project/quinnlab/dcl3nd/norfolk/data/_scratch_zarrs/" "/project/quinnlab/dcl3nd/norfolk/data/mrms_for_rainyday_hourly/" "/project/quinnlab/dcl3nd/norfolk/data/mrms_for_rainyday_daily/"

#%% testing on rivanna interactive job from python terminal
# module purge
# module load anaconda
# source activate mrms_processing
# in_date = 20210503

# f_in_nc = "/project/quinnlab/dcl3nd/norfolk/data/mrms_for_rainyday/" + "{}.nc".format(in_date)
# fl_out_zar = "/project/quinnlab/dcl3nd/norfolk/data/_scratch_zarrs/" + "{}_hourly.zarr".format(in_date)
# f_out_nc_hrly = "/project/quinnlab/dcl3nd/norfolk/data/mrms_for_rainyday_hourly/" + "{}.nc".format(in_date)
# f_out_nc_daily = "/project/quinnlab/dcl3nd/norfolk/data/mrms_for_rainyday_daily/" + "{}.nc".format(in_date)
#%% # inputs
in_date = str(sys.argv[1]) # YYYYMMDD
if "NULL" in in_date:
    sys.exit("Failed to create netcdf for {}. No netcdf file created likely because the SLURM_ARRAY_TASK_ID is 366 on a non-leap year. This is expected.".format(in_date))
f_in_nc = str(sys.argv[2]) + "{}.nc".format(in_date)
fl_out_zar = str(sys.argv[3]) + "{}_hourly.zarr".format(in_date)
f_out_nc_hrly = str(sys.argv[4]) + "{}.nc".format(in_date)
f_out_nc_daily = str(sys.argv[5]) + "{}.nc".format(in_date)

#%% load dataset
try:
    ds = xr.open_dataset(f_in_nc, chunks = {"time":chnk_sz})
except:
    files = glob(f_in_nc)
    if len(files)>0:
        print("ERROR: The file is present but the dataset failed to load. The matching file is:")
        print(files)
    else:
        print("A file does not exist for this day.")
        sys.exit()

#%% create dataset of hourly data
# resample to hourly
# (since data is already in mm/hr, a simple average will result in total mm of precipitation)
ds_hourly = ds.resample(time='1H').mean(skipna=True) # taking the mean preserves NA values

ds_hourly.rainrate.attrs["long_name"] = "Total Hourly Precipitation"
ds_hourly.rainrate.attrs["short_name"] = "tot_precip"
ds_hourly.rainrate.attrs["units"] = "mm"
ds_hourly.rainrate.attrs["description"] = "Radar total hourly precipitation"

ds_hourly.attrs = ds.attrs
del ds_hourly.attrs['rainrate_units']
ds_hourly.attrs["original_time_step"] = ds_hourly.attrs["time_step"] 
del ds_hourly.attrs["time_step"] 

#%% export to zarr
bm_time = time.time()
ds_hourly.to_zarr(fl_out_zar, mode="w")
# print("Created zarr: {}".format(time.time() - bm_time))

#%% load zarr and export to netcdf

ds_from_zarr = xr.open_zarr(store=fl_out_zar, chunks={'time':chnk_sz})
ds_from_zarr.to_netcdf(f_out_nc_hrly, encoding= {"rainrate":{"zlib":True}})
shutil.rmtree(fl_out_zar)
print("Created hourly netcdf: {}".format(time.time() - bm_time))



#%% load data
ds_hourly = xr.open_dataset(f_out_nc_hrly, chunks = {"time":chnk_sz})

#%% resample to daily timestep
# convert from mm/hr to mm/day
ds_hourly['rainrate'] = ds_hourly.rainrate * 24 # 24 hours per day

ds_daily = ds_hourly.resample(time='24H').mean(skipna=True)

ds_daily.rainrate.attrs["long_name"] = "Total Daily Precipitation"
ds_daily.rainrate.attrs["short_name"] = "tot_precip"
ds_daily.rainrate.attrs["units"] = "mm"
ds_daily.rainrate.attrs["description"] = "Radar total daily precipitation"

ds_daily.attrs = ds_hourly.attrs

#%% export to netcdf
bm_time = time.time()
ds_daily_loaded = ds_daily.load()
# ds_from_zarr = xr.open_zarr(store=fl_out_zar, chunks={'time':chnk_sz})
ds_daily_loaded.to_netcdf(f_out_nc_daily, encoding= {"rainrate":{"zlib":True}})
# print("Created daily netcdf: {}".format(time.time() - bm_time))

#%% finish
print("Script complete. Daily and hourly netcdfs created. Total runtime: {}".format(time.time() - start_time))