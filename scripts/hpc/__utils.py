#%% import libraries
import pandas as pd
import sys

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

# target chunk size per script
da_chnk_sz = "10000MB" # script seems to succeed when this is 1/4 of the memory allocated (I didn't try to push it)
db_chnk_sz = "5000MB"
dc_chnk_sz = "5000MB" 
ha_chnk_sz = "5000MB" 
hb_chnk_sz = "5000MB" 
i_chnk_sz = "10000MB"
i_chnk_sz_space = 90 # determined through trial and error measuring completion speed
#%% functions
def remove_vars(ds, coords_to_delete, attrs_to_delete):
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

def return_chunking_parameters(script_prefix):
    if script_prefix == "da":
        chnk_sz = da_chnk_sz
    elif script_prefix == "db":
        chnk_sz = db_chnk_sz
    elif script_prefix == "dc":
        chnk_sz = dc_chnk_sz
    elif script_prefix == "ha":
        chnk_sz = ha_chnk_sz
    elif script_prefix == "hb":
        chnk_sz = hb_chnk_sz
    elif script_prefix == "i":
        chnk_sz = i_chnk_sz
    else:
        sys.exit("The chunks size for this script has not been assigned!")

    return chnk_sz, size_of_float32, MB_per_bit, num_lats, num_lons









