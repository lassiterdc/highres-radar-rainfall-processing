#%% libraries and directories
import xarray as xr
from glob import glob
from tqdm import tqdm
import numpy as np
import dask
import pandas as pd
import shutil
import time
import geopandas as gpd
import rioxarray
import matplotlib.pyplot as plt
import sys
from __utils import return_chunking_parameters
from __utils import return_colorbar_percentile_for_plotting_gridded_precip_data
import __utils
import time
dask.config.set(**{'array.slicing.split_large_chunks': True})
bm_time = time.time()

coords_to_delete = ["step", "heightAboveSea", "valid_time"] # do not contain useful information
attrs_to_delete = ['source', 'problems'] # not valid for aggregated timestep

chnk_sz, size_of_float32, MB_per_bit, num_lats, num_lons = return_chunking_parameters("hb")
cbar_percentile, nearest_int_for_rounding = return_colorbar_percentile_for_plotting_gridded_precip_data()

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

#%% work 
f_in_nc_yearlyavg = "/project/quinnlab/dcl3nd/norfolk/highres-radar-rainfall-processing/data/mrms_nc_preciprate_yearly_singlefile.nc" # str(sys.argv[1])
fl_states = "/project/quinnlab/dcl3nd/norfolk/highres-radar-rainfall-processing/data/geospatial/States_shapefile.shp" #str(sys.argv[2])
f_shp_nexrad_boundary = "/project/quinnlab/dcl3nd/norfolk/highres-radar-rainfall-processing/arcpro/shapefiles/nexrad_boundary.shp" # str(sys.argv[3])
fldr_plots = "/project/quinnlab/dcl3nd/norfolk/highres-radar-rainfall-processing/plots/h_annual_statistics/" + "{}.png" # str(sys.argv[4]) + "{}.png"
#%% load input parameters
f_in_nc_yearlyavg = str(sys.argv[1])
fl_states = str(sys.argv[2])
f_shp_nexrad_boundary = str(sys.argv[3])
fldr_plots = str(sys.argv[4]) + "{}.png"

#%% load data
ds_yearly = xr.open_dataset(f_in_nc_yearlyavg)
bm_time = time.time()
ds_yearly = ds_yearly.load()

ds_yearly.rio.write_crs("epsg:4326", inplace=True)
ds_yearly.rio.set_spatial_dims("longitude", "latitude", inplace=True)

# convert latitude measurements to degrees east
ds_yearly["longitude"] = ds_yearly["longitude"] - 360


gdf_nexrad_boundary = gpd.read_file(f_shp_nexrad_boundary)
gdf_states = gpd.read_file(fl_states)
gdf_states = gdf_states[gdf_states.State_Name != "HAWAII"]
gdf_states = gdf_states[gdf_states.State_Name != "ALASKA"]

# clip to nexrad boundary
ds_yearly = ds_yearly.rio.clip(gdf_nexrad_boundary.geometry)

# compute range for colorbar
q_low = (1 - cbar_percentile) / 2
q_high = 1 -    q_low
highend_magnitude = abs(ds_yearly.rainrate.quantile(q_high).values)
lowend_magnitude = abs(ds_yearly.rainrate.quantile(q_low).values)
cbar_magnitude = max([highend_magnitude, lowend_magnitude])

cbar_magnitude = np.ceil(cbar_magnitude / nearest_int_for_rounding) * nearest_int_for_rounding
#%% plotting
#%% defining plot parameters
time_idxs = np.arange(len(ds_yearly.time))
ncols = int(np.ceil((width_to_height*len(time_idxs))**(0.5)/width_to_height)*width_to_height)
nrows = np.ceil(len(time_idxs) / ncols).astype(int)
#%% plotting rainfall totals
vmin = -1 * cbar_magnitude
vmax = cbar_magnitude

fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=[width, width/width_to_height], sharex = True, sharey = True)
time_idx = -1
for row in np.arange(nrows):
    for col in np.arange(ncols):
        time_idx += 1
        if time_idx > max(time_idxs):
            axes[row, col].axis('off')
            continue
        # plot data
        im = ds_yearly.rainrate.isel(time=time_idx).plot.pcolormesh(x="longitude", y="latitude", ax=axes[row, col], #col="time"
                                        cmap='jet', # col_wrap = ncols
                                        vmin = vmin, vmax = vmax, 
                                        # cbar_kwargs={"label":"This scale is colored according to the 2nd and 98th percentiles."},
                                        add_colorbar = False)
        # add shapefile plot
        gdf_states.plot(ax=axes[row, col], edgecolor='black', color="none", zorder=100)
        # add single colorbar
        if time_idx ==  max(time_idxs):
            cbar_ax = fig.add_axes([0.92, 0.13, 0.01, 0.7])
            fig.colorbar(im, cax=cbar_ax, pad=0.02, shrink=0.5, label="Annual Precipitation Total (mm) (colored based on {}th percnetile)".format(str(int(cbar_percentile*100))))
        # label graphs
        if col == 0 and row == 1:
            axes[row, col].set_ylabel("Latitude")
        else:
            axes[row, col].set_ylabel("")
        if row+1 == nrows and col == 2:
            axes[row, col].set_xlabel("Longitude")
        else:
            axes[row, col].set_xlabel("")
        axes[row, col].set_title(str(int(ds_yearly.time[time_idx].values)))
# ds_yearly.rainrate.plot.pcolormesh(x="longitude", y="latitude", col="time",
#                                    col_wrap = ncols, robust=True, figsize = [width, width/width_to_height], cmap='jet',
#                                    cbar_kwargs={"label":"This scale is colored according to the 2nd and 98th percentiles."})

plt.savefig(fldr_plots.format("all_years"), dpi=300)

#%% Temporary addition for 2023 ESE symposium for recruitment weekend
# ds_yearly.rainrate.plot.pcolormesh(x="longitude", y="latitude", col="time",
#                                    col_wrap = ncols, robust=True, figsize = [7.91, 7.91], cmap='jet',
#                                    cbar_kwargs={"label":"This scale is colored according to the 2nd and 98th percentiles."})
# plt.savefig(fldr_plots.format("all_years_ese_colloq_1"), dpi=300)

# ds_yearly.rainrate.plot.pcolormesh(x="longitude", y="latitude", col="time",
#                                    col_wrap = ncols, robust=True, figsize = [7.91, 6], cmap='jet',
#                                    cbar_kwargs={"label":"This scale is colored according to the 2nd and 98th percentiles."})
# plt.savefig(fldr_plots.format("all_years_ese_colloq_2"), dpi=300)



#%% plotting anomolies
# ds_yearly_rainrate = ds_yearly.rainrate.copy()
# # if quantized data is not being used, exclude these years
# idx_times = []
# if exclude_2013to2014 == True:
#     exclude_times = [2013, 2014]
#     for year in ds_yearly.time.values:
#         if year not in exclude_times:
#             idx_times.append(year)
#     ds_yearly_rainrate = ds_yearly.rainrate.sel({"time":idx_times})

# ds_yearly_rainrate_mean = ds_yearly_rainrate.mean(dim=["time"])

# ds_yearly_anomolies = ds_yearly.rainrate - ds_yearly_rainrate_mean

# ds_yearly_anomolies.plot.pcolormesh(x="longitude", y="latitude", col="time",
#                                    col_wrap = ncols, robust=True, figsize = [width, width/width_to_height], cmap='jet',
#                                    cbar_kwargs={"label":"This scale is colored according to the 2nd and 98th percentiles."})
# plt.savefig(fldr_plots.format("all_years_anomolies_rel_to_avg_excluding_2012-2014"), dpi=300)

print("Script runtime: {}".format(time.time() - bm_time))

