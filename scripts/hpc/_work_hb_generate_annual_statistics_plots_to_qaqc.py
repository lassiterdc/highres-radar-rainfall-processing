#%% code for testing on local machine
# environment: mrms_analysis
f_in_nc_dailyavg = 'D:/Dropbox/_GradSchool/_norfolk/highres-radar-rainfall-processing/data/mrms_daily_totals.nc'
f_in_nc_yearlyavg = 'D:/Dropbox/_GradSchool/_norfolk/highres-radar-rainfall-processing/data/mrms_yearly_totals.nc'
fl_states = "D:/Dropbox/_GradSchool/_norfolk/highres-radar-rainfall-processing/data/geospatial/States_shapefile.shp"
fldr_plots = "D:/Dropbox/_GradSchool/_norfolk/highres-radar-rainfall-processing/plots/h_annual_statistics/{}.png"

#%% code for testing in interactive job
# ijob -c 1 -A quinnlab_paid -p standard --time=1-00:00:00 --mem=33000
# module purge
# module load gcc openmpi eccodes anaconda
# source activate mrms_analysis
# norfolk
# python

# f_in_nc = "data/mrms_for_rainyday_daily_consolidated/*.nc"
# f_out_nc_dailyavg = 'data/mrms_daily_totals.nc'
# f_in_nc_yearlyavg = 'data/mrms_yearly_totals.nc'
# fl_out_zar = '_scratch_zarrs/h_yearly.zarr'
# fl_states = "data/shapefiles/States_shapefile.shp"
# fldr_plots = "plots/h_annual_statistics/{}.png"

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

dask.config.set(**{'array.slicing.split_large_chunks': True})

coords_to_delete = ["step", "heightAboveSea", "valid_time"] # do not contain useful information
attrs_to_delete = ['source', 'problems'] # not valid for aggregated timestep

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
#%% load input parameters
# f_in_ncs = str(sys.argv[1])
# shp_gages = str(sys.argv[2])
# f_out_csv = str(sys.argv[3])
# f_out_nc = str(sys.argv[4])

#%% functions
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

# #%% load data
# files = glob(f_in_nc)
# files.sort()
# lst_ds = []
# days = [] 

# for f in tqdm(files):
#     ds = xr.open_dataset(f, chunks={"latitude":chnk_lat, "longitude":chnk_lon})
#     ds = ds.sortby(["time"])
#     ds = remove_vars(ds, coords_to_delete, attrs_to_delete)
#     lst_ds.append(ds)
#     days.append(len(ds.time))

# ds_allyrs = xr.concat(lst_ds, dim="time", coords='minimal')
# #%% resample to year
# # https://pandas.pydata.org/pandas-docs/stable/user_guide/timeseries.html#time-date-components
# # ds_allyrs = ds_allyrs.sortby(["time"])
# ds_yearly = ds_allyrs.resample(time='Y').mean(skipna=True) 

# dtime = pd.DatetimeIndex(ds_yearly.time.values, freq='Y').year

# ds_yearly['time'] = dtime

# ds_yearly.rainrate.attrs["long_name"] = "Average Daily Precipitation"
# ds_yearly.rainrate.attrs["short_name"] = "avg_daily_precip"
# ds_yearly.rainrate.attrs["units"] = "mm/day"
# ds_yearly.rainrate.attrs["description"] = "Radar average daily precipitation"

# ds_yearly.attrs = ds.attrs
# # del ds_yearly.attrs['source']

# # del ds_yearly.attrs['problems']

# ds_yearly = ds_yearly.chunk(chunks={"longitude":chnk_lon, "latitude":chnk_lat, "time":1})
# # ds_yearly.attrs["original_time_step"] = ds.attrs["time_step"] 
# # del ds_yearly.attrs["time_step"] 

# #%% export
# bm_time = time.time()
# ds_yearly.to_zarr(fl_out_zar, mode="w")
# ds_from_zarr = xr.open_zarr(store=fl_out_zar, chunks={'time':chnk_sz})
# ds_from_zarr.to_netcdf(f_out_nc_dailyavg, encoding= {"rainrate":{"zlib":True}})
# shutil.rmtree(fl_out_zar)
# print("Created netcdf of annual daily averages: {}".format(time.time() - bm_time))


# #%% exporting yearly totals
# ds_yearly = xr.open_dataset(f_out_nc_dailyavg)
# bm_time = time.time()
# ds_yearly = ds_yearly.load()
# print("Loaded netcdf into memory: {}".format(time.time() - bm_time))

# #%% export annual totals
# # convert from mm/day to mm/year
# ds_yearly['rainrate'] = ds_yearly.rainrate * 365.25 # mm/day * days per year = mm/year

# #%% update attributes
# ds_yearly.rainrate.attrs["long_name"] = "Average Annual Precipitation"
# ds_yearly.rainrate.attrs["short_name"] = "avg_yearly_precip"
# ds_yearly.rainrate.attrs["units"] = "mm/year"
# ds_yearly.rainrate.attrs["description"] = "Radar average annual precipitation"

# #%% export
# bm_time = time.time()
# ds_yearly.to_netcdf(f_out_nc_yearlyavg, encoding= {"rainrate":{"zlib":True}})
# print("Created netcdf of annual totals: {}".format(time.time() - bm_time))


#%% load_dataset
ds_yearly = xr.open_dataset(f_in_nc_yearlyavg)
bm_time = time.time()
ds_yearly = ds_yearly.load()
print("Loaded netcdf into memory: {}".format(time.time() - bm_time))


######## plotting ###############
#%% defining plot parameters
width = 14
width_to_height = 1.5
tot_graphs = np.arange(len(ds_yearly.time))
ncols = int(np.ceil((width_to_height*len(tot_graphs))**(0.5)/width_to_height)*width_to_height)

#%% plotting rainfall totals
ds_yearly.rainrate.plot.pcolormesh(x="longitude", y="latitude", col="time",
                                   col_wrap = ncols, robust=True, figsize = [width, width/width_to_height], cmap='jet',
                                   cbar_kwargs={"label":"This scale is colored according to the 2nd and 98th percentiles."})
# plt.tight_layout()
plt.savefig(fldr_plots.format("all_years"), dpi=300)

#%% plotting anomolies
idx_times = []
exclude_times = [2012, 2013, 2014]
for year in ds_yearly.time.values:
    if year not in exclude_times:
        idx_times.append(year)

# exclude 2012-2015
ds_yearly_rainrate_subset = ds_yearly.rainrate.sel({"time":idx_times})
ds_yearly_rainrate_subset_mean = ds_yearly_rainrate_subset.mean(dim=["time"])

ds_yearly_anomolies = ds_yearly.rainrate - ds_yearly_rainrate_subset_mean

ds_yearly_anomolies.plot.pcolormesh(x="longitude", y="latitude", col="time",
                                   col_wrap = ncols, robust=True, figsize = [width, width/width_to_height], cmap='jet',
                                   cbar_kwargs={"label":"This scale is colored according to the 2nd and 98th percentiles."})
# plt.tight_layout()
plt.savefig(fldr_plots.format("all_years_anomolies_rel_to_avg_excluding_2012-2014"), dpi=300)


#%% virginia
# https://stackoverflow.com/questions/51398563/python-mask-netcdf-data-using-shapefile
from rasterio import features
from affine import Affine

def transform_from_latlon(lat, lon):
    """ input 1D array of lat / lon and output an Affine transformation
    """
    lat = np.asarray(lat)
    lon = np.asarray(lon)
    trans = Affine.translation(lon[0], lat[0])
    scale = Affine.scale(lon[1] - lon[0], lat[1] - lat[0])
    return trans * scale

def rasterize(shapes, coords, latitude='latitude', longitude='longitude',
              fill=np.nan, **kwargs):
    """Rasterize a list of (geometry, fill_value) tuples onto the given
    xray coordinates. This only works for 1d latitude and longitude
    arrays.

    usage:
    -----
    1. read shapefile to geopandas.GeoDataFrame
          `states = gpd.read_file(shp_dir+shp_file)`
    2. encode the different shapefiles that capture those lat-lons as different
        numbers i.e. 0.0, 1.0 ... and otherwise np.nan
          `shapes = (zip(states.geometry, range(len(states))))`
    3. Assign this to a new coord in your original xarray.DataArray
          `ds['states'] = rasterize(shapes, ds.coords, longitude='X', latitude='Y')`

    arguments:
    ---------
    : **kwargs (dict): passed to `rasterio.rasterize` function

    attrs:
    -----
    :transform (affine.Affine): how to translate from latlon to ...?
    :raster (numpy.ndarray): use rasterio.features.rasterize fill the values
      outside the .shp file with np.nan
    :spatial_coords (dict): dictionary of {"X":xr.DataArray, "Y":xr.DataArray()}
      with "X", "Y" as keys, and xr.DataArray as values

    returns:
    -------
    :(xr.DataArray): DataArray with `values` of nan for points outside shapefile
      and coords `Y` = latitude, 'X' = longitude.


    """
    transform = transform_from_latlon(coords[latitude], coords[longitude])
    out_shape = (len(coords[latitude]), len(coords[longitude]))
    raster = features.rasterize(shapes, out_shape=out_shape,
                                fill=fill, transform=transform,
                                dtype=float, **kwargs)
    spatial_coords = {latitude: coords[latitude], longitude: coords[longitude]}
    return xr.DataArray(raster, coords=spatial_coords, dims=(latitude, longitude))

def add_shape_coord_from_data_array(xr_da, shp_path, coord_name):
    """ Create a new coord for the xr_da indicating whether or not it 
         is inside the shapefile

        Creates a new coord - "coord_name" which will have integer values
         used to subset xr_da for plotting / analysis/

        Usage:
        -----
        precip_da = add_shape_coord_from_data_array(precip_da, "awash.shp", "awash")
        awash_da = precip_da.where(precip_da.awash==0, other=np.nan) 
    """
    # 1. read in shapefile
    shp_gpd = gpd.read_file(shp_path)

    # 2. create a list of tuples (shapely.geometry, id)
    #    this allows for many different polygons within a .shp file (e.g. States of US)
    shapes = [(shape, n) for n, shape in enumerate(shp_gpd.geometry)]

ds_yearly = add_shape_coord_from_data_array(ds_yearly, fl_states, "virginia")
da_va = ds_yearly.where(ds_yearly.virginia==0, other=np.nan)
da_va.mean(dim="time").plot()