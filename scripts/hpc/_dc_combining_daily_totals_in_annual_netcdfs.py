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
import pathlib
from __utils import return_chunking_parameters


chnk_sz = return_chunking_parameters("dc")[0]
remove_daily_files = False # set to true if you want to remove all the individual netcdfs for each day (set to false while still developing script)
#%% # inputs
in_date = str(sys.argv[1]) # YYYYMMDD
if "NULL" in in_date:
    sys.exit("Failed to create netcdf for {}. No netcdf file created likely because the SLURM_ARRAY_TASK_ID is 366 on a non-leap year. This is expected.".format(in_date))
f_in_nc = str(sys.argv[2]) + "{}*.nc".format(in_date)
fl_out_zar = str(sys.argv[3]) + "{}_yearly.zarr".format(in_date)
f_out_nc_daily = str(sys.argv[4]) + "{}.nc".format(in_date)

#%% load dataset
try:
    ds = xr.open_mfdataset(f_in_nc, concat_dim="time", combine='nested',
                     combine_attrs="drop_conflicts", chunks={"time":1},
                     coords="minimal")
except:
    files = glob(f_in_nc)
    # if len(files)>0:
    #     print("ERROR: The file is present but the dataset failed to load. The matching file is:")
    #     print(files)
    # else:
    #     print("Files does not exist for this day.")
    sys.exit("Failed to open dataset for {}. Filename: {}".format(in_date, files))

ds = xr.open_mfdataset(f_in_nc, concat_dim="time", combine='nested',
                     combine_attrs="drop_conflicts", chunks={"time":1},
                     coords="minimal")

# print("Loaded dataset. Creating zarr file...")
#%% export to zarr
bm_time = time.time()
ds.to_zarr(fl_out_zar, mode="w")

# load zarr and export to netcdf
ds_from_zarr = xr.open_zarr(store=fl_out_zar, chunks={'time':chnk_sz})
# print("Time to create zarr: {}".format(time.time() - bm_time))
ds_from_zarr.to_netcdf(f_out_nc_daily, encoding= {"rainrate":{"zlib":True}})

# delete zarr file
shutil.rmtree(fl_out_zar)
# print("Created year-long netcdf of daily totals by first exporting to zarr then to netcdf: {}".format(time.time() - bm_time))
#%% remove individual netcdfs
if remove_daily_files == True:
    files = glob(f_in_nc)
    for f in files:
        pth = pathlib.Path(f)
        pth.unlink()
#%% finish
print("Daily and hourly netcdfs created for {}. Total runtime: {}".format(in_date, time.time() - start_time))