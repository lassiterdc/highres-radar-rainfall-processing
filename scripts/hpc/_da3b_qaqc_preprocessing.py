#%% Import libraries
import shutil
import pandas as pd
import sys
from glob import glob
import xarray as xr
import time
from __utils import *
import os

start_time = time.time()

#%% work
fldr_out_nc = "/scratch/dcl3nd/highres-radar-rainfall-processing/data/mrms_nc_preciprate_fullres_dailyfiles_constant_tstep/"
fldr_csvs = "/scratch/dcl3nd/highres-radar-rainfall-processing/data/_scratch/csv/"
fl_nc_qaqc_out_all = fldr_out_nc + "_qaqc_of_resampled_data.nc"
grouping_index = 0

#%% end work
fldr_out_nc = str(sys.argv[1]) # ${assar_dirs[out_fullres_dailyfiles_consolidated]} # "/scratch/dcl3nd/highres-radar-rainfall-processing/out_fullres_dailyfiles_consolidated/"
grouping_index = int(sys.argv[2])-1 # subtracting 1 since the slurm task id is 1-indexed
fl_nc_qaqc_out_all = fldr_out_nc + "_qaqc_of_resampled_data.nc"

#%% create more aggregated summary statistic netcdfs
bm_time = time.time()
import dask
dask.config.set(**{'array.slicing.split_large_chunks': False})
ds_qaqc_all = xr.open_dataset(fl_nc_qaqc_out_all, engine="h5netcdf", chunks = {'date': -1})

# create list of grouping variables
ds_qaqc_all.coords["year_month"] = ds_qaqc_all['date'].dt.strftime('%Y-%m')
groupings = ['date.year', 'date.month', "year_month"]
f_pattern = fldr_out_nc + "_qaqc_of_resampled_data_{}_{}.nc"

def write_netcdf(ds_out, f_out):
    bm_time_intermediate1 = time.time() 
    # write zarr file 
    f_zarr = f_out.split(".nc")[0] + ".zarr"
    ds_out.to_zarr(f_zarr, mode="w")
    time_elapsed_min = round((time.time() - bm_time_intermediate1) / 60, 2)
    print("{} minutes to export export zarr".format(time_elapsed_min))
    # load zarr file
    ds_out = xr.open_dataset(f_zarr, engine = "zarr")
    bm_time_intermediate = time.time() 
    # write netcdf file
    ds_out.to_netcdf(f_out, engine = "h5netcdf")
    time_elapsed_min = round((time.time() - bm_time_intermediate) / 60, 2)
    print("{} additional minutes to export export netcdf".format(time_elapsed_min))
    shutil.rmtree(f_zarr)
    time_elapsed_min = round((time.time() - bm_time_intermediate1) / 60, 2)
    print("{} minutes total for writing {}".format(time_elapsed_min, f_out))
    return

group_var =  groupings[grouping_index]

# min
f_out = f_pattern.format(group_var, "min")
ds_out = ds_qaqc_all.groupby(group_var).min(dim = "date")#.load()
write_netcdf(ds_out, f_out)
# max
f_out = f_pattern.format(group_var, "max")
ds_out = ds_qaqc_all.groupby(group_var).max("date")
write_netcdf(ds_out, f_out)
# mean
f_out = f_pattern.format(group_var, "mean")
ds_out = ds_qaqc_all.groupby(group_var).mean("date")
write_netcdf(ds_out, f_out)
#quantiles
f_out = f_pattern.format(group_var, "quants")
ds_out = ds_qaqc_all.groupby(group_var).quantile(q = [0.1,0.5,0.9], dim = "date", method = "hazen")
write_netcdf(ds_out, f_out)

time_elapsed_min = round((time.time() - bm_time) / 60, 2)
print("Exported consolidated qaqc netcdf files. Time to export: {}.".format(time_elapsed_min))

    