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
from __utils import *
#%% parameters
chnk_sz = db_chnk_sz
d_perf = {}
#%% testing

in_date = "20020719" #str(sys.argv[1]) # YYYYMMDD
f_in_nc = "/scratch/dcl3nd/highres-radar-rainfall-processing/data/mrms_nc_preciprate_fullres_dailyfiles_constant_tstep/" + "{}.nc".format(in_date)
fl_out_zar = "/scratch/dcl3nd/highres-radar-rainfall-processing/data/_scratch/zarrs/" + "{}_hourly.zarr".format(in_date)
f_out_nc_hrly = "/scratch/dcl3nd/highres-radar-rainfall-processing/data/mrms_nc_preciprate_hourly_dailyfiles/"+ "{}.nc".format(in_date)
f_out_nc_daily = "/scratch/dcl3nd/highres-radar-rainfall-processing/data/mrms_nc_preciprate_daily_dailyfiles/" + "{}.nc".format(in_date)
f_out_csv = "/scratch/dcl3nd/highres-radar-rainfall-processing/data/_scratch/csv/" + "db_consolidating_tseps_{}.csv".format(in_date)

#%% # inputs
in_date = str(sys.argv[1]) # YYYYMMDD
f_in_nc = str(sys.argv[2]) + "{}.nc".format(in_date)
fl_out_zar = str(sys.argv[3]) + "{}_hourly.zarr".format(in_date)
f_out_nc_hrly = str(sys.argv[4]) + "{}.nc".format(in_date)
f_out_nc_daily = str(sys.argv[5]) + "{}.nc".format(in_date)
f_out_csv = str(sys.argv[6]) + "db_consolidating_tseps_{}.csv".format(in_date)

d_perf["f_in_nc"] = f_in_nc
success = False
if ("NULL" not in in_date) and ("qaqc" not in f_in_nc.lower()): # if the date is valid
    d_perf["date"] = int(in_date)
    #%% load dataset
    success = True
    try:
        # load dataset keeping only the rainrate variable
        ds = xr.Dataset()
        ds["rainrate"] = xr.open_dataset(f_in_nc, chunks = {"time":chnk_sz}).rainrate
        # delete unused coordinates 
        ds = ds.reset_coords(drop=True)
    except Exception as e:
        success = False
        d_perf["error_running_xr.open_dataset()"] = e

    #%% create dataset of hourly data
    # resample to hourly
    # (since data is already in mm/hr, a simple average will result in total mm of precipitation)
    if success: # if the dataset opened succesfully, resample
        ds_hourly = ds.resample(time='1H').mean(skipna=True) # taking the mean preserves NA values
        ds_hourly.rainrate.attrs["long_name"] = "Total Hourly Precipitation"
        ds_hourly.rainrate.attrs["short_name"] = "tot_precip"
        ds_hourly.rainrate.attrs["units"] = "mm"
        ds_hourly.rainrate.attrs["description"] = "Radar total hourly precipitation"
        ds_hourly.attrs = ds.attrs
        try:
            del ds_hourly.attrs['rainrate_units']
        except:
            pass
        try:
            ds_hourly.attrs["original_time_step"] = ds_hourly.attrs["time_step"] 
        except:
            pass
        try:
            del ds_hourly.attrs["time_step"]
        except:
            pass 
        #%% export to zarr
        bm_time = time.time()
        success = True
        try:
            ds_hourly.to_zarr(fl_out_zar, mode="w")
        except Exception as e:
            success = False
            d_perf["error_running_to_zarr)"] = e
        #%% load zarr and export to netcdf
        if success: # if the dataset exported to zarr
            success = True
            try:
                ds_from_zarr = xr.open_zarr(store=fl_out_zar, chunks={'time':chnk_sz})
                ds_from_zarr.to_netcdf(f_out_nc_hrly, encoding= {"rainrate":{"zlib":True}})
                shutil.rmtree(fl_out_zar)
            except Exception as e:
                success = False
                d_perf["error_running_to_netcdf_on_hrly"] = e
            # print("Created hourly netcdf: {}".format(time.time() - bm_time))
            #%% load data
            if success: # if the hourly dataset was exported to a netcdf
                success = True
                try:
                    ds_hourly = xr.open_dataset(f_out_nc_hrly, chunks = {"time":chnk_sz})
                    # resample to daily timestep
                    # convert from mm/hr to mm/day
                    ds_hourly['rainrate'] = ds_hourly.rainrate * 24 # 24 hours per day

                    ds_daily = ds_hourly.resample(time='24H').mean(skipna=True)

                    ds_daily.rainrate.attrs["long_name"] = "Total Daily Precipitation"
                    ds_daily.rainrate.attrs["short_name"] = "tot_precip"
                    ds_daily.rainrate.attrs["units"] = "mm"
                    ds_daily.rainrate.attrs["description"] = "Radar total daily precipitation"

                    ds_daily.attrs = ds_hourly.attrs
                    # export to netcdf
                    bm_time = time.time()
                    ds_daily_loaded = ds_daily.load()
                    ds_daily_loaded.to_netcdf(f_out_nc_daily, encoding= {"rainrate":{"zlib":True}})
                    d_perf["total_runtime"] = time.time() - start_time
                except Exception as e:
                    success = False
                    d_perf["error_running_to_netcdf_on_daily"] = e
    d_perf["success"] = success

# print(d_perf)
df = pd.DataFrame(d_perf, index = [0])
df.to_csv(f_out_csv)