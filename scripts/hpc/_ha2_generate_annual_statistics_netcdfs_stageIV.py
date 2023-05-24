#%% libraries and directories
import xarray as xr
from glob import glob
from tqdm import tqdm
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

chnk_sz, size_of_float32, MB_per_bit, num_lats, num_lons = return_chunking_parameters("ha2")

days_in_year = 365
total_size_MB = days_in_year * num_lats * num_lons * size_of_float32 * MB_per_bit
target_chunks_size_MB = int(chnk_sz.split("MB")[0])
num_chunks = total_size_MB / target_chunks_size_MB
chnks_per_dim = np.sqrt(num_chunks)
chnk_lat = int(round(num_lats / chnks_per_dim))
chnk_lon = int(round(num_lons / chnks_per_dim))

#%% load input parameters
# work
f_in_nc = "/project/quinnlab/dcl3nd/norfolk/stormy/data/climate/StageIV_rainfall/" + "2003/*20030501*.nc" # str(sys.argv[1])
f_out_nc_yearlyavg = "/project/quinnlab/dcl3nd/norfolk/highres-radar-rainfall-processing/data/stageiv_nc_preciprate_yearly_singlefile.nc" # str(sys.argv[2])
fl_out_zar = "/project/quinnlab/dcl3nd/norfolk/highres-radar-rainfall-processing/data/_scratch/zarrs/" + 'ha2_yearly.zarr' # str(sys.argv[3])
fl_states = "/project/quinnlab/dcl3nd/norfolk/highres-radar-rainfall-processing/data/geospatial/States_shapefile.shp" # str(sys.argv[4]) 
fldr_plots = "/project/quinnlab/dcl3nd/norfolk/highres-radar-rainfall-processing/plots/h_annual_statistics/" + "{}.png" # str(sys.argv[4]) + "{}.png"
# end work

f_in_nc = str(sys.argv[1]) + "*/*.nc" # "/project/quinnlab/dcl3nd/norfolk/stormy/data/climate/StageIV_rainfall/"
f_out_nc_yearlyavg = str(sys.argv[2]) # "data/stageiv_nc_preciprate_yearly_singlefile.nc"
fl_out_zar = str(sys.argv[3]) + 'ha2_yearly.zarr' # "data/_scratch/zarrs/"
fl_states = str(sys.argv[4]) # "data/geospatial/States_shapefile.shp"

#%% load data
files = glob(f_in_nc)
files.sort()
lst_ds = []
days = []

# strt = time.time()
# for f in tqdm(files):
# # for f in files:
#     ds = xr.open_dataset(f, chunks={"latitude":chnk_lat, "longitude":chnk_lon}, engine='h5netcdf')
#     # ds = ds.sortby(["time"])
#     ds = remove_vars(ds)
#     lst_ds.append(ds)
#     days.append(len(ds.time))
#     end = time.time()
# print("Created netcdf of annual averages. Script runtime: {}".format(end - strt))
# ds_allyrs = xr.concat(lst_ds, dim="time", coords='minimal')
# ds_allyrs = ds_allyrs.sortby(["time"])


# strt = time.time()
ds_allyrs = xr.open_mfdataset(files, chunks={"latitude":chnk_lat, "longitude":chnk_lon}, engine = "h5netcdf")
# convert from preceding to following time interval
ds_allyrs["time"] = ds_allyrs.time - pd.Timedelta(1, "h")

## convert any negative values to nan
ds_allyrs = ds_allyrs.where(ds_allyrs["rainrate"]>=0)

## convert any values above 9000 to nan
ds_allyrs = ds_allyrs.where(ds_allyrs["rainrate"]<9000)

# print("Created netcdf of annual averages. Script runtime: {}".format(time.time() - strt))

#%% work - test plotting
# ds_allyrs.rainrate.isel(time=[0,1]).plot.pcolormesh(x='outlon', y='outlat', col="time",
#                                    robust=False, cmap='jet')
#                                 #    vmin = 0, vmax = 1800)
#                                    #cbar_kwargs={"label":"This scale is colored according to the 2nd and 98th percentiles."})
# plt.savefig(fldr_plots.format("stage_iv_test1"), dpi=300)
#%%


#%% resample to year
# https://pandas.pydata.org/pandas-docs/stable/user_guide/timeseries.html#time-date-components
ds_yearly = ds_allyrs.resample(time='Y').mean(skipna=True) 

#%% work - test plotting
# ds_yearly.rainrate.plot.pcolormesh(x='outlon', y='outlat', col="time",
#                                    robust=False, cmap='jet')
#                                 #    vmin = 0, vmax = 1800)
#                                    #cbar_kwargs={"label":"This scale is colored according to the 2nd and 98th percentiles."})
# plt.savefig(fldr_plots.format("stage_iv_test1"), dpi=300)
#%%

dtime = pd.DatetimeIndex(ds_yearly.time.values, freq='Y').year

ds_yearly['time'] = dtime

ds_yearly.rainrate.attrs["long_name"] = "Average Yearly Precipitation"
ds_yearly.rainrate.attrs["short_name"] = "avg_yearly_precip"
ds_yearly.rainrate.attrs["units"] = "mm/year"
ds_yearly.rainrate.attrs["description"] = "Radar average yearly precipitation"

# ds_yearly.attrs = ds.attrs

ds_yearly = ds_yearly.chunk(chunks={"outlat":chnk_lon, "outlon":chnk_lat, "time":1})

ind_time = []
lst_ds = []

# convert from mm/day to mm/year
for yr in ds_yearly.time.values:
    if yr % 4 == 0: # it is a leap year
        hours = 366 * 24
    else:
        hours = 365 * 24
    ds = ds_yearly.rainrate.loc[dict(time=yr)] * hours # mm/hr * hours/year = mm/year
    lst_ds.append(ds)
da_allyrs_mmperyear = xr.concat(lst_ds, dim="time", coords='minimal')
ds_yearly["rainrate"] = da_allyrs_mmperyear
#%% work
# ds_yearly.to_netcdf(f_out_nc_yearlyavg, encoding= {"rainrate":{"zlib":True}})
# print("Created netcdf of annual averages. Script runtime: {}".format(time.time() - bm_time))
#%% export
ds_yearly.to_zarr(fl_out_zar, mode="w")
ds_from_zarr = xr.open_zarr(store=fl_out_zar, chunks={'time':chnk_sz})
ds_from_zarr.to_netcdf(f_out_nc_yearlyavg, encoding= {"rainrate":{"zlib":True}})

shutil.rmtree(fl_out_zar)
print("Created netcdf of annual averages. Script runtime: {}".format(time.time() - bm_time))