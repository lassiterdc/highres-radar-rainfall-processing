#%% Import libraries and set filepaths
import xarray as xr
import numpy as np
import pandas as pd
import time
import dask
dask.config.set(**{'array.slicing.split_large_chunks': False})

f_in_ncs = "D:/mrms_processing/data/mrms_for_rainyday_subset_norfolk_netcdfs/*.nc"
f_out_nc = "D:/mrms_processing/out_netcdfs/a1_mrms_at_norfolk.nc"
f_out_nc_24h = "D:/mrms_processing/out_netcdfs/a1_mrms_at_norfolk_24h.nc"
f_out_nc_1M = "D:/mrms_processing/out_netcdfs/a1_mrms_at_norfolk_1M.nc"
#%% load netcdf data
mrms = xr.open_mfdataset(f_in_ncs,  concat_dim = "time",
            combine = "nested", engine = 'netcdf4', coords='minimal')

attrs_to_drop = ["source", "warnings", "missing_timesteps_filled_with_NA_values", 
"coordinates_adjusted_to_conform_to_modern_grid", "previous_latitude_values_on_grid_centers",
"previous_longitude_values_on_grid_centers", 
"mean_latitude_shift_in_translation_to_modern_grid (most recent minus previous, so positive means the most recent grid was more north)",
"mean_longitude_shift_in_translation_to_modern_grid (most recent minus previous, so positive means the most recent grid was more east)",
]

vars_to_drop = ['step', 'heightAboveSea', 'valid_time']

for at in attrs_to_drop:
    try:
        del mrms.attrs[at]
    except:
        continue

for var in vars_to_drop:
    try:
        del mrms[var]
    except:
        continue

#%% summarize to other timesteps
def add_attributes(new_ds, tstep, desc, units, old_ds = mrms):
    new_ds.rainrate.attrs["long_name"] = "{} precipitation {} ({})".format(tstep, desc, units)
    new_ds.rainrate.attrs["short_name"] = "{}_precip_{}".format(tstep, units)
    new_ds.rainrate.attrs["units"] = units
    new_ds.rainrate.attrs["description"] = "radar {} precipitation".format(tstep)

    new_ds.attrs = old_ds.attrs
    return new_ds

def convert_mm_per_hr_to_mm_per_month(ds):
    hrs_in_month = pd.to_datetime(ds.time.values).days_in_month * 24
    ds["rainrate"] = ds.rainrate * hrs_in_month # mm/hr * hrs
    return ds

def convert_mm_per_hr_to_mm_per_day(ds):
    hrs_in_day = 24
    ds["rainrate"] = ds.rainrate * hrs_in_day # mm/hr * hrs
    return ds
#%%
mrms_sorted = mrms.sortby("time")

# monthly
print("Creating monthly netcdf timeseries...")
bm_mnth = time.time()
mrms_1m_mm_per_hr = mrms_sorted.resample(time="1M").mean()
mrms_1m = mrms_1m_mm_per_hr.groupby("time").map(convert_mm_per_hr_to_mm_per_month)
mrms_1m = add_attributes(mrms_1m, tstep = "monthly", desc="totals", units="mm")
mrms_rsmpld_loaded = mrms_1m.load()
mrms_rsmpld_loaded.to_netcdf(f_out_nc_1M)
bm_mnth = time.time() - bm_mnth

# daily
print("Creating daily netcdf timeseries...")
bm_day = time.time()
mrms_1d_mm_per_hour = mrms_sorted.resample(time="24H").mean()
mrms_1d = mrms_1d_mm_per_hour.groupby("time").map(convert_mm_per_hr_to_mm_per_day)
mrms_1d = add_attributes(mrms_1d, tstep = "daily", desc="totals", units="mm")
mrms_1d_loaded = mrms_1d.load()
mrms_1d_loaded.to_netcdf(f_out_nc_24h)
bm_day = time.time() - bm_day

# export full dataset
print("Exporting netdf with all data...")
bm_full = time.time()
mrms_loaded = mrms.load()
mrms_loaded.to_netcdf(f_out_nc)
mrms_loaded.close()
bm_full = time.time() - bm_full

print("Finished exporting netcdfs. BMs: {} for month, {} for day, {} for all".format(bm_mnth, bm_day, bm_full))