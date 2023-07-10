#%% Import libraries
# import time
# start_time = time.time()
# import shutil
import xarray as xr
# import cfgrib
# from glob import glob
# import numpy as np
# from scipy import stats
import pandas as pd
import sys
# from tqdm import tqdm
import dask
# import os
dask.config.set(**{'array.slicing.split_large_chunks': False}) # to silence warnings of loading large slice into memory
dask.config.set(scheduler='synchronous') # this forces single threaded computations
# from pathlib import Path
# import pathlib
# import __utils
# from __utils import remove_vars
# from __utils import return_corner_coords
from __utils import return_chunking_parameters
# from __utils import return_target_tstep

target_tstep = 5

chnk_sz = "10MB"

#%% work
# fldr_in_nc_day = "/project/quinnlab/dcl3nd/norfolk/highres-radar-rainfall-processing/data/mrms_nc_preciprate_fullres_dailyfiles/"
# fldr_out_nc = "/scratch/dcl3nd/highres-radar-rainfall-processing/mrms_nc_preciprate_fullres_dailyfiles_{}min/".format(target_tstep)
# in_date = "20200605"
#%% end work

# folders (with proceeding fwd slash)
in_date = str(sys.argv[1]) # YYYYMMDD
fldr_in_nc_day = str(sys.argv[2]) # ${assar_dirs[out_fullres_dailyfiles]} # "/project/quinnlab/dcl3nd/norfolk/highres-radar-rainfall-processing/data/mrms_nc_preciprate_fullres_dailyfiles/"
fldr_out_nc = str(sys.argv[3]) # ${assar_dirs[out_fullres_dailyfiles_5min]} # "/scratch/dcl3nd/highres-radar-rainfall-processing/mrms_nc_preciprate_fullres_dailyfiles_5min/"

if "NULL" in in_date:
    sys.exit("Failed to create netcdf for {}. No netcdf file created likely because the SLURM_ARRAY_TASK_ID is 366 on a non-leap year. This is expected.".format(in_date))

#%% netcdf 
fl_in_nc = fldr_in_nc_day +"{}.nc".format(in_date)
fl_out_nc = fldr_out_nc +"{}.nc".format(in_date)

# print(fl_in_nc)

ds = xr.open_dataset(fl_in_nc, chunks = dict(latitude = chnk_sz))

# tstep = ds.attrs["time_step"]
tstep_min = pd.to_timedelta(ds.attrs["time_step"]).total_seconds() / 60

num_tsteps = ds.coords["time"].shape[0]

duration_h = num_tsteps * tstep_min / 60

# verify the full day has coverage
if duration_h != 24:
    sys.exit("Failed to create netcdf for {}. Duration does not add up to 24 hours.".format(in_date))

#%% testing 
# da = xr.DataArray(np.arange(100), coords=[pd.date_range("2020-01-01", periods = 100, freq = "2T")], dims = "time")

# t_idx_1min = pd.date_range("2020-01-01", periods = 200, freq='1min')

# da_1min = da.reindex(dict(time = t_idx_1min)).ffill(dim="time")

# da_5min = da_1min.resample(time = "5T").mean()

#%% resampling
# da = xr.DataArray(np.arange(100), coords=[pd.date_range("2020-01-01", periods = 100, freq = "2T")], dims = "time")

t_idx_1min = pd.date_range(ds.time.values[0], periods = 24*60*60, freq='1min')

ds_1min = ds.reindex(dict(time = t_idx_1min)).ffill(dim="time")

da_target = ds_1min.resample(time = "{}T".format(target_tstep)).mean()

da_target = da_target.chunk(chunks = dict(latitude = chnk_sz))  

da_target = da_target.unify_chunks()

da_target.to_netcdf(fl_out_nc, encoding= {"rainrate":{"zlib":True}})