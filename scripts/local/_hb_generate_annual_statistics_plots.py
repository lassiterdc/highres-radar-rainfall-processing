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
from __filepaths import return_chunking_parameters
from __filepaths import return_colorbar_percentile_for_plotting_gridded_precip_data
from __filepaths import *
import __filepaths
import time
dask.config.set(**{'array.slicing.split_large_chunks': True})
# Set matplotlib parameters
plt.rcParams['text.usetex'] = True
# Set the axes title font size
plt.rc('axes', titlesize=16)
# Set the font size of the figure title
plt.rc('figure', titlesize=16)
# Set the axes labels font size
plt.rc('axes', labelsize=14)
# Set the font size for x tick labels
plt.rc('xtick', labelsize=12)
# Set the font size for y tick labels
plt.rc('ytick', labelsize=12)
# Set the legend font size
plt.rc('legend', fontsize=16)

#%% load and preprocess data for plotting
bm_time = time.time()

coords_to_delete = ["step", "heightAboveSea", "valid_time"] # do not contain useful information
attrs_to_delete = ['source', 'problems'] # not valid for aggregated timestep

chnk_sz, size_of_float32, MB_per_bit, num_lats, num_lons = return_chunking_parameters("hb")
cbar_percentile, nearest_int_for_rounding = return_colorbar_percentile_for_plotting_gridded_precip_data()

# use_quantized_data = __utils.use_quantized_data
exclude_2013to2014 = __filepaths.exclude_2013to2014_from_mean_for_anamolies_plot
# plotting parameters
width = __filepaths.plt_hb_width
width_to_height = __filepaths.plt_hb_width / __filepaths.plt_hb_height

days_in_year = 365
total_size_MB = days_in_year * num_lats * num_lons * size_of_float32 * MB_per_bit
target_chunks_size_MB = int(chnk_sz.split("MB")[0])
num_chunks = total_size_MB / target_chunks_size_MB
chnks_per_dim = np.sqrt(num_chunks)
chnk_lat = int(round(num_lats / chnks_per_dim))
chnk_lon = int(round(num_lons / chnks_per_dim))

# load data
ds_yearly = xr.open_dataset(f_nc_yearly)
bm_time = time.time()

ds_yearly.rio.write_crs("epsg:4326", inplace=True)
ds_yearly.rio.set_spatial_dims("longitude", "latitude", inplace=True)

# convert latitude measurements to degrees east
ds_yearly["longitude"] = ds_yearly["longitude"] - 360


gdf_nexrad_boundary = gpd.read_file(f_shp_nexrad_boundary)
gdf_states = gpd.read_file(f_shp_states)
gdf_states = gdf_states[gdf_states.State_Name != "HAWAII"]
gdf_states = gdf_states[gdf_states.State_Name != "ALASKA"]

# defining plot parameters
time_idxs = np.arange(len(ds_yearly.time))
ncols = 5
nrows = np.ceil(len(time_idxs) / ncols).astype(int)

df_dim_counts = pd.DataFrame((np.unique(np.concatenate((np.diff(ds_yearly.latitude.values), np.diff(ds_yearly.longitude.values))), return_counts=True))).T
df_dim_counts.columns = ['dim_spacing', 'count']
grid_resolution = np.around(df_dim_counts.iloc[df_dim_counts['count'].idxmax(), 0], 2) # round to nearest hundredth of a degree
coarsen_factor = int(coarsen_to / grid_resolution) # aggregation of lat and long to achieve desired grid resolution

coarsen = True
if coarsen == True:
    ds_yearly_plt = ds_yearly.coarsen(dim=dict(latitude = coarsen_factor, longitude = coarsen_factor)).mean()

else:
    ds_yearly_plt = ds_yearly

ds_yearly_plt = ds_yearly_plt.rio.clip(gdf_nexrad_boundary.geometry)
ds_yearly_plt = ds_yearly_plt.load()

# convert from mm to inches
ds_yearly_plt['rainrate'] = ds_yearly_plt.rainrate / mm_per_in

# compute range for colorbar
# q_low = (1 - cbar_percentile) / 2
q_high = cbar_percentile
highend_magnitude = abs(ds_yearly_plt.rainrate.quantile(q_high).values)
# lowend_magnitude = abs(ds_yearly.rainrate.quantile(q_low).values)
cbar_magnitude = highend_magnitude

cbar_magnitude = np.ceil(cbar_magnitude / nearest_int_for_rounding) * nearest_int_for_rounding
#%% plotting rainfall totals
vmin = 0
vmax = cbar_magnitude

cbar_str = str(int(cbar_percentile*100))

fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=[width*0.9, width/width_to_height*0.75],
                          sharex = True, sharey = True, gridspec_kw = {'wspace':0.05, 'hspace':0.1},
                          dpi=600)
time_idx = -1
for row in np.arange(nrows):
    for col in np.arange(ncols):
        time_idx += 1
        if time_idx > max(time_idxs):
            axes[row, col].axis('off')
            continue
        # plot data
        im = ds_yearly_plt.rainrate.isel(time=time_idx).plot.pcolormesh(x="longitude", y="latitude", ax=axes[row, col], #col="time"
                                        cmap='jet_r', # col_wrap = ncols
                                        vmin = vmin, vmax = vmax, 
                                        # cbar_kwargs={"label":"This scale is colored according to the 2nd and 98th percentiles."},
                                        add_colorbar = False)
        # # add shapefile plot
        gdf_states.plot(ax=axes[row, col], edgecolor='black', color="none", zorder=100)
        # add single colorbar
        
        if time_idx ==  max(time_idxs):
            # im = ds_yearly_plt.rainrate.isel(time=time_idx).plot.pcolormesh(x="longitude", y="latitude", ax=axes[row, col], #col="time"
            #                                 cmap='jet', # col_wrap = ncols
            #                                 vmin = vmin, vmax = vmax, 
            #                                 # cbar_kwargs={"label":"This scale is colored according to the 2nd and 98th percentiles."},
            #                                 add_colorbar = False)
            cbar_ax = fig.add_axes([0.92, 0.13, 0.01, 0.7])
            fig.colorbar(im, cax=cbar_ax, pad=0.02, shrink=0.5, label="MRMS Annual Totals (in)")
        # label graphs
        # if col == 0 and row+1 == nrows :
        #     axes[row, col].set_xlabel("Longitude")
        #     axes[row, col].set_ylabel("Latitude")
        # else:
        #     axes[row, col].set_xlabel("")
        #     axes[row, col].set_ylabel("")
        axes[row, col].set_xlabel("")
        axes[row, col].set_ylabel("")
        if row+1 == nrows:
            labs = axes[row, col].get_xticklabels()
            ticks = axes[row, col].get_xticks()
            count = -1
            for lab in labs:
                count += 1
                old_text = lab.get_text()
                val = int(ticks[count])
                new_text = str(-1 * val) + "$^{\circ}$W"
                lab.set_text(new_text)

            axes[row, col].set_xticklabels(labs)

        if col == 0:
            labs = axes[row, col].get_yticklabels()
            ticks = axes[row, col].get_yticks()
            count = -1
            for lab in labs:
                count += 1
                old_text = lab.get_text()
                val = int(ticks[count])
                new_text = str(val) + "$^{\circ}$N"
                lab.set_text(new_text)

            axes[row, col].set_yticklabels(labs)


        axes[row, col].set_title(str(int(ds_yearly.time[time_idx].values)))
print('saving figure....')
plt.savefig(fldr_out_plots_h + "mrms_annual_totals_cbar{}_{}degree_spacing.png".format(cbar_str, coarsen_to))