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
from __utils import return_chunking_parameters
import __utils
import time
dask.config.set(**{'array.slicing.split_large_chunks': True})
bm_time = time.time()

coords_to_delete = ["step", "heightAboveSea", "valid_time"] # do not contain useful information
attrs_to_delete = ['source', 'problems'] # not valid for aggregated timestep

chnk_sz, size_of_float32, MB_per_bit, num_lats, num_lons = return_chunking_parameters("hb")

# use_quantized_data = __utils.use_quantized_data
exclude_2013to2014 = __utils.exclude_2013to2014_from_mean_for_anamolies_plot
# plotting parameters
width = __utils.plt_hb_width
width_to_height = __utils.plt_hb_width / __utils.plt_hb_height

days_in_year = 365
total_size_MB = days_in_year * num_lats * num_lons * size_of_float32 * MB_per_bit
target_chunks_size_MB = int(chnk_sz.split("MB")[0])
num_chunks = total_size_MB / target_chunks_size_MB
chnks_per_dim = np.sqrt(num_chunks)
chnk_lat = int(round(num_lats / chnks_per_dim))
chnk_lon = int(round(num_lons / chnks_per_dim))

#%% load input parameters
f_in_nc_yearlyavg = str(sys.argv[1])
fl_states = str(sys.argv[2])
fldr_plots = str(sys.argv[3]) + "{}.svg"

#%% load_dataset
ds_yearly = xr.open_dataset(f_in_nc_yearlyavg)
bm_time = time.time()
ds_yearly = ds_yearly.load()
# print("Loaded netcdf into memory: {}".format(time.time() - bm_time))
#%% plotting
#%% defining plot parameters
tot_graphs = np.arange(len(ds_yearly.time))
ncols = int(np.ceil((width_to_height*len(tot_graphs))**(0.5)/width_to_height)*width_to_height)

#%% plotting rainfall totals
ds_yearly.rainrate.plot.pcolormesh(x="longitude", y="latitude", col="time",
                                   col_wrap = ncols, robust=True, figsize = [width, width/width_to_height], cmap='jet',
                                   cbar_kwargs={"label":"This scale is colored according to the 2nd and 98th percentiles."})
plt.savefig(fldr_plots.format("all_years"), dpi=300)

#%% Temporary addition for 2023 ESE symposium for recruitment weekend
ds_yearly.rainrate.plot.pcolormesh(x="longitude", y="latitude", col="time",
                                   col_wrap = ncols, robust=True, figsize = [7.91, 7.91], cmap='jet',
                                   cbar_kwargs={"label":"This scale is colored according to the 2nd and 98th percentiles."})
plt.savefig(fldr_plots.format("all_years_ese_colloq_1"), dpi=300)

ds_yearly.rainrate.plot.pcolormesh(x="longitude", y="latitude", col="time",
                                   col_wrap = ncols, robust=True, figsize = [7.91, 6], cmap='jet',
                                   cbar_kwargs={"label":"This scale is colored according to the 2nd and 98th percentiles."})
plt.savefig(fldr_plots.format("all_years_ese_colloq_2"), dpi=300)



#%% plotting anomolies
ds_yearly_rainrate = ds_yearly.rainrate.copy()
# if quantized data is not being used, exclude these years
idx_times = []
if exclude_2013to2014 == True:
    exclude_times = [2013, 2014]
    for year in ds_yearly.time.values:
        if year not in exclude_times:
            idx_times.append(year)
    ds_yearly_rainrate = ds_yearly.rainrate.sel({"time":idx_times})

ds_yearly_rainrate_mean = ds_yearly_rainrate.mean(dim=["time"])

ds_yearly_anomolies = ds_yearly.rainrate - ds_yearly_rainrate_mean

ds_yearly_anomolies.plot.pcolormesh(x="longitude", y="latitude", col="time",
                                   col_wrap = ncols, robust=True, figsize = [width, width/width_to_height], cmap='jet',
                                   cbar_kwargs={"label":"This scale is colored according to the 2nd and 98th percentiles."})
plt.savefig(fldr_plots.format("all_years_anomolies_rel_to_avg_excluding_2012-2014"), dpi=300)

print("Script runtime: {}".format(time.time() - bm_time))

