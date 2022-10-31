#%% libraries and directories
import xarray as xr
from glob import glob
# from tqdm import tqdm
import numpy as np
import dask
import pandas as pd
import shutil
import time
import matplotlib.pyplot as plt
import sys
import __utils
from __utils import remove_vars
from __utils import return_chunking_parameters
dask.config.set(**{'array.slicing.split_large_chunks': True})
bm_time = time.time()

include_2012_2013_and_2014 = __utils.use_quantized_data

chnk_sz, size_of_float32, MB_per_bit, num_lats, num_lons = return_chunking_parameters("ha")

days_in_year = 365
total_size_MB = days_in_year * num_lats * num_lons * size_of_float32 * MB_per_bit
target_chunks_size_MB = int(chnk_sz.split("MB")[0])
num_chunks = total_size_MB / target_chunks_size_MB
chnks_per_dim = np.sqrt(num_chunks)
chnk_lat = int(round(num_lats / chnks_per_dim))
chnk_lon = int(round(num_lons / chnks_per_dim))

#%% load input parameters
f_in_nc = str(sys.argv[1]) + "*.nc"
f_out_nc_yearlyavg = str(sys.argv[2])
fl_out_zar = str(sys.argv[3]) + 'h_yearly.zarr'
fl_states = str(sys.argv[4])

#%% load data
files = glob(f_in_nc)
files.sort()
lst_ds = []
days = []

for f in files:
    ds = xr.open_dataset(f, chunks={"latitude":chnk_lat, "longitude":chnk_lon})
    ds = ds.sortby(["time"])
    ds = remove_vars(ds)
    lst_ds.append(ds)
    days.append(len(ds.time))

ds_allyrs = xr.concat(lst_ds, dim="time", coords='minimal')
#%% resample to year
# https://pandas.pydata.org/pandas-docs/stable/user_guide/timeseries.html#time-date-components
ds_yearly = ds_allyrs.resample(time='Y').mean(skipna=True) 

dtime = pd.DatetimeIndex(ds_yearly.time.values, freq='Y').year

ds_yearly['time'] = dtime

ds_yearly.rainrate.attrs["long_name"] = "Average Yearly Precipitation"
ds_yearly.rainrate.attrs["short_name"] = "avg_yearly_precip"
ds_yearly.rainrate.attrs["units"] = "mm/year"
ds_yearly.rainrate.attrs["description"] = "Radar average yearly precipitation"

ds_yearly.attrs = ds.attrs

ds_yearly = ds_yearly.chunk(chunks={"longitude":chnk_lon, "latitude":chnk_lat, "time":1})

ind_time = []
if include_2012_2013_and_2014 == False:
    for t in dtime:
        if t not in [2012, 2013, 2014]:
            ind_time.append(t)
    ds_yearly = ds_yearly.sel({"time":ind_time})

# convert from mm/day to mm/year
for yr in ds_yearly.time.values:
    if yr % 4 == 0: # it is a leap year
        days = 366
    else:
        days = 365
    ds_yearly.rainrate.loc[dict(time=yr)] = ds_yearly.rainrate.loc[dict(time=yr)] * days # mm/day * days/year = mm/year

#%% export

ds_yearly.to_zarr(fl_out_zar, mode="w")
ds_from_zarr = xr.open_zarr(store=fl_out_zar, chunks={'time':chnk_sz})
ds_from_zarr.to_netcdf(f_out_nc_yearlyavg, encoding= {"rainrate":{"zlib":True}})


#%% remove zarr file
shutil.rmtree(fl_out_zar)
print("Created netcdf of annual averages. Script runtime: {}".format(time.time() - bm_time))