#%% import libraries
import pandas as pd
import sys


#%% filepaths
# fldr_repo = "/project/quinnlab/dcl3nd/norfolk/highres-radar-rainfall-processing/"
fldr_repo = "/scratch/dcl3nd/highres-radar-rainfall-processing/"
fldr_nc_fullres_daily = fldr_repo + "data/mrms_nc_preciprate_fullres_dailyfiles/"
fldr_nc_fullres_daily_constant_tstep = fldr_repo + "data/mrms_nc_preciprate_fullres_dailyfiles_constant_tstep/"
fldr_scratch_zarr = fldr_repo + "data/_scratch/zarrs/"
fldr_scratch_csv = fldr_repo + "data/_scratch/csv/"
f_shp_sst_transom = fldr_repo + "/stormy/stochastic_storm_transposition/norfolk/transposition_domain/norfolk_trans_dom_4326.shp"
fldr_nc_stageiv = "/scratch/dcl3nd/stormy/data/climate/StageIV_rainfall/"
#%% user options
# use data from 8/2013 - 2014?
# this will likely stay 'False' because I found that
# these data truncated heavy rainfall causing obvious
# underestimation of annual rainfall totals
use_quantized_data = False

# gage data information
gage_id_attribute_in_shapefile = "MONITORING"

# plotting options
## script hb
plt_hb_width = 14
plt_hb_height = plt_hb_width / 1.5
# even 2013-2014 data is included in analysis, it is preferable not to include it in the computation of mean annual rainfall
exclude_2013to2014_from_mean_for_anamolies_plot = True
#%% hard coding
# associated function: remove_vars
coords_to_delete = ["step", "heightAboveSea", "valid_time"] # do not contain useful information
attrs_to_delete = ['source', 'problems'] # not valid for aggregated timestep

# These coordiantes represent the latitude and longitude coordinates
# of the northwest and southeast most gridcells in the dataset
# provided by Dr. Jian Zhang, the lead on the 2001-2011 MRMS reanalysis
# dataset (https://doi.org/10.25638/EDC.PRECIP.0001)
s_nw_corner_lat = 54.995 # degrees north
s_nw_corner_lon = 129.995 # degrees west
s_se_corner_lat = 20.005 # degrees north
s_se_corner_lon = 60.005 # degrees west

# chunking parameters
size_of_float32 = 32 # bits
MB_per_bit = 1.25e-7
num_lats = 3500
num_lons = 7000

lst_quants = [0.1,0.5,0.9] # quantiles of rainfall correction factors for each day of consolidated rainyday ready mrms data
crxn_upper_bound = 20 # upper bound of correction factor
crxn_lower_bound = 0.01 # lower bound of correction factor
# target timestep
target_tstep = 5 # minutes

# target chunk size per script
da_chnk_sz = "10000MB" # script seems to succeed when this is 1/4 of the memory allocated (I didn't try to push it)
db_chnk_sz = "5000MB"
dc_chnk_sz = "5000MB" 
ha_chnk_sz = "5000MB" 
hb_chnk_sz = "5000MB" 
i_chnk_sz = "10000MB"
i_chnk_sz_space = 90 # determined through trial and error measuring completion speed

# percentile for colorbar for plotting mrms data
cbar_percentile = 0.98
nearest_int_for_rounding = 50
#%% functions
def process_dans_stageiv(ds_st4):
    import xarray as xr
    import numpy as np
    ds_st4['outlat'] = ds_st4.latitude.values
    ds_st4['outlon'] = ds_st4.longitude.values+360
    ds_st4 = ds_st4.drop_vars("latitude")
    ds_st4 = ds_st4.drop_vars("longitude")
    ds_st4 = ds_st4.drop_vars("infilled")
    ds_st4 = ds_st4.rename({"outlat":"latitude", "outlon":"longitude"})
    # replace negative values with np.nan
    ds_st4 = xr.where(ds_st4>0, ds_st4, np.nan) # where condition is true, keep ds_stageiv; else fill with np.nan
    return ds_st4

def clip_ds_to_another_ds(ds_to_clip, ds_target, lat_varname="latitude", lon_varname="longitude"):
    # assuming ds_target is rectangular
    ## define coordinate boundaries based on ds_target
    max_lat = ds_target[lat_varname].values.max()
    min_lat = ds_target[lat_varname].values.min()
    max_lon = ds_target[lon_varname].values.max()
    min_lon = ds_target[lon_varname].values.min()
    ## subset ds_to_clip based on boundaries defined above
    ds_clipped = ds_to_clip.where((ds_to_clip[lat_varname] >= min_lat) & (ds_to_clip[lat_varname] <= max_lat) & (ds_to_clip[lon_varname] >= min_lon) & (ds_to_clip[lon_varname] <= max_lon), drop = True)
    return ds_clipped

def spatial_resampling(xds_to_resample, xds_target, lat_varname, lon_varname, missingfillval = 0):
    from rasterio.enums import Resampling
    import xarray as xr
    # load the dataset with the target resolution
    # rename dimensions to x and y
    xds_target = xds_target.rename({lat_varname:"y", lon_varname:"x"})
    # load the dataset to be modified
    xds_to_resample = xds_to_resample.rename({lat_varname:"y", lon_varname:"x"})
    # assign coordinate system
    xds_target.rio.write_crs("epsg:4326", inplace=True)
    xds_to_resample.rio.write_crs("epsg:4326", inplace=True)
    # set spatial dimensions
    xds_target.rio.set_spatial_dims("x", "y", inplace=True)
    xds_to_resample.rio.set_spatial_dims("x", "y", inplace=True)
    # resample
    ## https://corteva.github.io/rioxarray/stable/rioxarray.html#rioxarray.raster_dataset.RasterDataset.reproject_match 
    ## (https://rasterio.readthedocs.io/en/stable/api/rasterio.enums.html#rasterio.enums.Resampling)
    xds_to_resampled = xds_to_resample.rio.reproject_match(xds_target, resampling = Resampling.average)
    # rename back to original dimension names
    xds_to_resampled = xds_to_resampled.rename({"y":lat_varname, "x":lon_varname})
    # fill missing values with prespecified val (this should just corresponds to areas where one dataset has pieces outside the other)
    xds_to_resampled = xr.where(xds_to_resampled>=3.403e+37, x = missingfillval, y = xds_to_resampled)
    return xds_to_resampled


# def return_target_tstep():
#     return target_tstep

def remove_vars(ds, coords_to_delete=coords_to_delete, attrs_to_delete=attrs_to_delete):
    '''
    This removes coordinates and attributes that are not useful for the combined MRMS product. These are
    hardcoded in the __utils.py file
    '''

    for crd in coords_to_delete:
        try:
            del ds.coords[crd]
        except: 
            continue
    for att in attrs_to_delete:
        try:
            del ds.attrs[att]
        except:
            continue
    return ds

def __convert_degrees_west_to_degrees_east(x):
    coord_conv = (360 - x)
    return coord_conv

def return_corner_coords():
    """Returns 4 hardcoded variables: the latitude and longitude of the northwest corner of the grid 
    and the latitude and longitude of the southeast corner of the grid, respecitvely.
    Latitude units are degrees north and longitude units are in degrees east"""
    s_nw_corner_lon_deg_east = __convert_degrees_west_to_degrees_east(s_nw_corner_lon)
    s_se_corner_lon_deg_east = __convert_degrees_west_to_degrees_east(s_se_corner_lon)
    return s_nw_corner_lat, s_nw_corner_lon_deg_east, s_se_corner_lat, s_se_corner_lon_deg_east

# def return_chunking_parameters(script_prefix, num_lats=num_lats, num_lons=num_lons):
#     if script_prefix == "da":
#         chnk_sz = da_chnk_sz
#     elif script_prefix == "db":
#         chnk_sz = db_chnk_sz
#     elif script_prefix == "dc":
#         chnk_sz = dc_chnk_sz
#     elif script_prefix == "ha":
#         chnk_sz = ha_chnk_sz
#     elif script_prefix == "ha2":
#         chnk_sz = ha_chnk_sz
#         num_lats = 900
#         num_lons = 2100
#     elif script_prefix == "hb":
#         chnk_sz = hb_chnk_sz
#     elif script_prefix == "i":
#         chnk_sz = i_chnk_sz
    # else:
    #     sys.exit("The chunks size for this script has not been assigned!")

    # return chnk_sz, size_of_float32, MB_per_bit, num_lats, num_lons

# def return_colorbar_percentile_for_plotting_gridded_precip_data():
#     return cbar_percentile, nearest_int_for_rounding








