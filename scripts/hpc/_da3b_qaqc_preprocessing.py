#%% Import libraries
import shutil
import pandas as pd
import sys
from glob import glob
import xarray as xr
import time
from __utils import *
import os
import flox
import flox.xarray

start_time = time.time()

#%% work
fldr_out_nc = "/scratch/dcl3nd/highres-radar-rainfall-processing/data/mrms_nc_preciprate_fullres_dailyfiles_constant_tstep/"
fldr_csvs = "/scratch/dcl3nd/highres-radar-rainfall-processing/data/_scratch/csv/"
fl_nc_qaqc_out_all = fldr_out_nc + "_qaqc_of_resampled_data.nc"
grouping_index = 2
# testing_quantile_only = False
#%% end work
fldr_out_nc = str(sys.argv[1]) # ${assar_dirs[out_fullres_dailyfiles_consolidated]} # "/scratch/dcl3nd/highres-radar-rainfall-processing/out_fullres_dailyfiles_consolidated/"
grouping_index = int(sys.argv[2])-1 # subtracting 1 since the slurm task id is 1-indexed
fl_nc_qaqc_out_all = fldr_out_nc + "_qaqc_of_resampled_data.nc"

#%% create more aggregated summary statistic netcdfs
bm_time = time.time()
import dask
dask.config.set(**{'array.slicing.split_large_chunks': False})
ds_qaqc_all = xr.open_dataset(fl_nc_qaqc_out_all, engine="h5netcdf", chunks = dict(latitude = 'auto', longitude = 'auto', date = -1))

# create list of grouping variables
ds_qaqc_all.coords["year_month"] = ds_qaqc_all['date'].dt.strftime('%Y-%m')
groupings = ['date.year', 'date.month', "year_month"]
f_pattern = fldr_out_nc + "_qaqc_of_resampled_data_{}_{}.nc"

def write_netcdf(ds_out, f_out, write_zarr_first = False):
    bm_time_intermediate1 = time.time()
    if write_zarr_first:
        # write zarr file 
        f_zarr = f_out.split(".nc")[0] + ".zarr"
        ds_out.to_zarr(f_zarr, mode="w")
        time_elapsed_min = round((time.time() - bm_time_intermediate1) / 60, 2)
        print("{} minutes to export export zarr".format(time_elapsed_min))
        # load zarr file
        ds_out = xr.open_dataset(f_zarr, engine = "zarr")
    bm_time_intermediate = time.time() 
    # write netcdf file
    d_encoding = {}
    for da_name in ds_out.data_vars:
        d_encoding[da_name] = {"zlib":True}
    ds_out.to_netcdf(f_out, encoding=d_encoding, engine = "h5netcdf")
    time_elapsed_min = round((time.time() - bm_time_intermediate) / 60, 2)
    if write_zarr_first:
        print("{} additional minutes to export export netcdf".format(time_elapsed_min))
        shutil.rmtree(f_zarr)
    time_elapsed_min = round((time.time() - bm_time_intermediate1) / 60, 2)
    print("{} minutes total for writing {}".format(time_elapsed_min, f_out))
    return

# rechunk for blockwise groupby operation
group_var =  groupings[grouping_index]
ds_qaqc_all = ds_qaqc_all.unify_chunks()
ds_qaqc_all = flox.rechunk_for_blockwise(ds_qaqc_all, axis = 'date', labels = ds_qaqc_all[group_var].values)

# if not testing_quantile_only:
# min
f_out = f_pattern.format(group_var, "min")
ds_out = ds_qaqc_all.groupby(group_var).min(dim = "date", method="blockwise", engine="flox")#.load()
write_netcdf(ds_out, f_out)
# max
f_out = f_pattern.format(group_var, "max")
ds_out = ds_qaqc_all.groupby(group_var).max("date", method="blockwise", engine="flox")
write_netcdf(ds_out, f_out)
# mean
f_out = f_pattern.format(group_var, "mean")
ds_out = ds_qaqc_all.groupby(group_var).mean("date", method="blockwise", engine="flox")
write_netcdf(ds_out, f_out)
#quantiles
f_out = f_pattern.format(group_var, "quants")
# ds = ds_rain_qaqc.groupby('date.year').min(dim = "date", method="blockwise", engine="flox")
# ds_out = ds_qaqc_all.groupby(group_var).quantile(q = [0.1,0.5,0.9], dim = "date", engine="flox", method="blockwise")
finalize_kwargs = dict(q=[0.1,0.5,0.9])
ds_out = flox.xarray.xarray_reduce(ds_qaqc_all, ds_qaqc_all[group_var], func="quantile",
                            method = "blockwise", engine = "flox", keep_attrs=True, **finalize_kwargs)
write_netcdf(ds_out, f_out)

time_elapsed_min = round((time.time() - bm_time) / 60, 2)
print("Exported consolidated qaqc netcdf files. Time to export: {}.".format(time_elapsed_min))

    
# %%