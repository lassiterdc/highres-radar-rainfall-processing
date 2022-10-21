#%% code for testing on local machine
# environment: mrms_analysis
# f_in_nc = "D:/mrms_processing/data/mrms_for_rainyday_daily_consolidated/*.nc"
# f_out_nc_dailyavg = 'D:/Dropbox/_GradSchool/_norfolk/highres-radar-rainfall-processing/data/mrms_daily_totals.nc'
# f_out_nc_yearlyavg = 'D:/Dropbox/_GradSchool/_norfolk/highres-radar-rainfall-processing/data/mrms_yearly_totals.nc'
# fl_out_zar = 'D:/Dropbox/_GradSchool/_norfolk/highres-radar-rainfall-processing/_temp/h_yearly.zarr'
# fl_states = "D:/Dropbox/_GradSchool/_norfolk/highres-radar-rainfall-processing/data/geospatial/States_shapefile.shp"
# fldr_plots = "D:/Dropbox/_GradSchool/_norfolk/highres-radar-rainfall-processing/plots/h_annual_statistics/{}.png"

#%% code for testing in interactive job
# ijob -c 1 -A quinnlab_paid -p standard --time=1-00:00:00 --mem=100000
# module purge
# module load gcc openmpi eccodes anaconda
# source activate mrms_analysis
# norfolk
# python

# f_in_nc = "data/mrms_for_rainyday_daily_consolidated/*.nc"
# f_out_nc_dailyavg = 'data/mrms_daily_totals.nc'
# f_out_nc_yearlyavg = 'data/mrms_yearly_totals.nc'
# fl_out_zar = '_scratch_zarrs/h_yearly.zarr'
# fl_states = "data/shapefiles/States_shapefile.shp"
# fldr_plots = "plots/h_annual_statistics/" + "{}.png"

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
from __utils import remove_vars

dask.config.set(**{'array.slicing.split_large_chunks': True})

# coords_to_delete = ["step", "heightAboveSea", "valid_time"] # do not contain useful information
# attrs_to_delete = ['source', 'problems'] # not valid for aggregated timestep

exclude_2012_2013_and_2014 = True

chnk_sz = "5000MB"
size_of_float32 = 32 # bits
MB_per_bit = 1.25e-7
num_lats = 3500
num_lons = 7000
tsteps = 365
total_size_MB = tsteps * num_lats * num_lons * size_of_float32 * MB_per_bit
target_chunks_size_MB = int(chnk_sz.split("MB")[0])
num_chunks = total_size_MB / target_chunks_size_MB
chnks_per_dim = np.sqrt(num_chunks)
chnk_lat = int(round(num_lats / chnks_per_dim))
chnk_lon = int(round(num_lons / chnks_per_dim))

#%% load input parameters
f_in_nc = str(sys.argv[1]) + "*.nc"
# f_out_nc_dailyavg = str(sys.argv[2])
f_out_nc_yearlyavg = str(sys.argv[2])
fl_out_zar = str(sys.argv[3]) + 'h_yearly.zarr'
fl_states = str(sys.argv[4])
# fldr_plots = str(sys.argv[5]) + "{}.png"

#%% functions (moved to __utils.py)
# def remove_vars(ds, coords_to_delete, attrs_to_delete):
#     for crd in coords_to_delete:
#         try:
#             del ds.coords[crd]
#         except: 
#             continue
#     for att in attrs_to_delete:
#         try:
#             del ds.attrs[att]
#         except:
#             continue
#     return ds

#%% load data
files = glob(f_in_nc)
files.sort()
lst_ds = []
days = []


for f in tqdm(files):
    ds = xr.open_dataset(f, chunks={"latitude":chnk_lat, "longitude":chnk_lon})
    ds = ds.sortby(["time"])
    ds = remove_vars(ds)
    lst_ds.append(ds)
    days.append(len(ds.time))

ds_allyrs = xr.concat(lst_ds, dim="time", coords='minimal')
#%% resample to year
# https://pandas.pydata.org/pandas-docs/stable/user_guide/timeseries.html#time-date-components
# ds_allyrs = ds_allyrs.sortby(["time"])
ds_yearly = ds_allyrs.resample(time='Y').mean(skipna=True) 

dtime = pd.DatetimeIndex(ds_yearly.time.values, freq='Y').year

ds_yearly['time'] = dtime

ds_yearly.rainrate.attrs["long_name"] = "Average Yearly Precipitation"
ds_yearly.rainrate.attrs["short_name"] = "avg_yearly_precip"
ds_yearly.rainrate.attrs["units"] = "mm/year"
ds_yearly.rainrate.attrs["description"] = "Radar average yearly precipitation"

ds_yearly.attrs = ds.attrs
# del ds_yearly.attrs['source']

# del ds_yearly.attrs['problems']

ds_yearly = ds_yearly.chunk(chunks={"longitude":chnk_lon, "latitude":chnk_lat, "time":1})
# ds_yearly.attrs["original_time_step"] = ds.attrs["time_step"] 
# del ds_yearly.attrs["time_step"] 

ind_time = []
if exclude_2012_2013_and_2014 == True:
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
bm_time = time.time()
ds_yearly.to_zarr(fl_out_zar, mode="w")
ds_from_zarr = xr.open_zarr(store=fl_out_zar, chunks={'time':chnk_sz})
ds_from_zarr.to_netcdf(f_out_nc_yearlyavg, encoding= {"rainrate":{"zlib":True}})
print("Created netcdf of annual averages: {}".format(time.time() - bm_time))

#%% remove zarr file
shutil.rmtree(fl_out_zar)