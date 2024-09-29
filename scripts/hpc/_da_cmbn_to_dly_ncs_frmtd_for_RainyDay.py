#%% Import libraries
import time
start_time = time.time()
import shutil
import xarray as xr
import cfgrib
from glob import glob
import numpy as np
from scipy import stats
import pandas as pd
import sys
from tqdm import tqdm
import dask
import os
dask.config.set(**{'array.slicing.split_large_chunks': False}) # to silence warnings of loading large slice into memory
# dask.config.set(scheduler='synchronous') # this forces single threaded computations
from pathlib import Path
import pathlib
import __utils
from __utils import remove_vars
from __utils import return_corner_coords
# from __utils import return_chunking_parameters
# importing hard coded variables
overwrite_previous = True

s_nw_corner_lat, s_nw_corner_lon, s_se_corner_lat, s_se_corner_lon = return_corner_coords()

# chnk_sz, size_of_float32, MB_per_bit, num_lats, num_lons = return_chunking_parameters("da")

use_quantized_data = __utils.use_quantized_data

## for testing
use_subset_of_files_for_testing = False  # this toggles the use of just the first 5 files of the dataset to be processed
if use_subset_of_files_for_testing:
    subset_size = 100

#%% work
# in_date = "20150619"
# fldr_mesonet_grib = "/project/quinnlab/dcl3nd/norfolk/highres-radar-rainfall-processing/data/raw_data/raw_data/mrms_grib_mesonet/" + "*{}*.grib2".format(in_date)
# fldr_nssl_grib = "/project/quinnlab/dcl3nd/norfolk/highres-radar-rainfall-processing/data/raw_data/raw_data/mrms_grib_nssl/" + "*{}*.grib2".format(in_date)
# fldr_mesonet_nc_frm_png = "/project/quinnlab/dcl3nd/norfolk/highres-radar-rainfall-processing/data/raw_data/mrms_nc_quant/" + "*{}*.nc".format(in_date)
# fldr_out_zar_day = "/project/quinnlab/dcl3nd/norfolk/highres-radar-rainfall-processing/data/_scratch/zarrs/"
# fldr_out_tmp_grib=" /project/quinnlab/dcl3nd/norfolk/highres-radar-rainfall-processing/data/_scratch/gribs/"
# fldr_out_nc_day = "/project/quinnlab/dcl3nd/norfolk/highres-radar-rainfall-processing/data/mrms_nc_preciprate_fullres_dailyfiles/"
# fldr_mesonet_grib_all = "/project/quinnlab/dcl3nd/norfolk/highres-radar-rainfall-processing/data/raw_data/raw_data/mrms_grib_mesonet/" + "*.grib2"
#%% end work

#%% # inputs
in_date = str(sys.argv[1]) # YYYYMMDD

# fail if this is for the 366th day of a non-leap year
if "NULL" in in_date:
    sys.exit("Failed to create netcdf for {}. No netcdf file created likely because the SLURM_ARRAY_TASK_ID is 366 on a non-leap year. This is expected.".format(in_date))

# folders (with proceeding fwd slash)
fldr_mesonet_grib = str(sys.argv[2]) + "*{}*.grib2".format(in_date)
fldr_nssl_grib = str(sys.argv[3]) + "*{}*.grib2".format(in_date)
fldr_mesonet_nc_frm_png = str(sys.argv[4]) + "*{}*.nc".format(in_date)
fldr_out_zar_day = str(sys.argv[5])
fldr_out_tmp_grib = str(sys.argv[6])
fldr_out_nc_day = str(sys.argv[7])
fldr_mesonet_grib_all = str(sys.argv[2]) + "*.grib2"

print(f"Processing MRMS data for date {in_date}")

#%%
fl_out_nc = fldr_out_nc_day +"{}.nc".format(in_date)
if use_subset_of_files_for_testing == True:
    fl_out_nc = fldr_out_nc_day +"{}_subset.nc".format(in_date)
# if the file exists and the script is set to NOT overwrite previous:
if os.path.isfile(fl_out_nc) and (not overwrite_previous):
    process_date = False
    print(f"File already exists. Not overwriting because overwrite_previous={overwrite_previous}. {fl_out_nc}")
    sys.exit("Stopping script...")
elif os.path.isfile(fl_out_nc) and (overwrite_previous):
    process_date = True
    print(f"File already exists but is being overwritten because overwrite_previous={overwrite_previous}. {fl_out_nc}")
else:
    process_date = True


# extract filename of most recent  from the Iowa State Environmental Mesonet
files_mesonet_grib_all = glob(fldr_mesonet_grib_all)
files_mesonet_grib_all.sort()
f_most_recent_grib = files_mesonet_grib_all[-1]

# extract most recent latitude and longitude
# print("Most recent file downloaded: {}".format(f_most_recent_grib))
with xr.open_dataset(f_most_recent_grib, engine='cfgrib') as ds:
    ds = ds.sortby(["latitude", "longitude"])
    v_lats_most_recent = ds["latitude"].astype(np.float32)
    v_lons_most_recent = ds["longitude"].astype(np.float32)
#%% files
files = glob(fldr_mesonet_grib) + glob(fldr_nssl_grib)

if use_quantized_data == True:
    files = files + glob(fldr_mesonet_nc_frm_png)

if use_subset_of_files_for_testing == True:
    files = files[0:subset_size]
    print("Warning: this script is running with a subset of the data to speedily inspect results. If this is undesired, set use_subset_of_files_for_testing to False.")

files.sort()
#%% function to extract timestep from filename
def extract_file_timestep(fname):
    fname = fname.split('/')[-1]
    fname = fname.split(".")
    ftype = fname.pop(-1)
    fname = ''.join(fname)
    str_tstep = fname.split("_")[-1]
    if ftype == "nc":
        date_format = '%Y%m%d%H%M'
    if ftype == "grib2":
        date_format = '%Y%m%d-%H%M%S'
    tstep = pd.to_datetime(str_tstep, format=date_format)
    return tstep
#%% Loading and formatting data
bm_time = time.time()

# load datasets and check to ensure that the latitude and longitude grids line up
## specify dictionary of warnings
dic_warnings = {"grid_is_unchanged_between_timesteps":True,
                "standard_timestep":True,
                "data_complete_ie_no_missing_timesteps":True,
                "data_conforms_to_standard_grid_dimensions":True,
                "notes_on_missing_values":"""Raw .grib2 data from the Iowa State Mesonet from 2015 onward uses -1 to 
flag missing values and -3 to flag 'fillvalues' which are areas outside of the radar domain. Based on analyzing
a small subset of raw data from the MRMS Reanalysis Dataset and then netcdfs dervied from quantized PNGs from the
Iowa State Mesonet, it appears that missing and fillvalues are assigned a value of 'nan'. Just to be safe,
the script for generating daily netcdfs from these data assigns a value of 'nan' to all negative values,
meaning a 'nan' can represent either missing data OR fillvalues."""}

lst_problems = []
lst_ds = []
# i = -1
if len(files) == 0: # if there's no data, stop script
    sys.exit("There appears to be no data for {} because the length of the list of files is 0.".format(in_date))

# to define chunks
tsteps = 720 # assuming 2 minute time intervals which will generate smaller chunks
# total_size_MB = tsteps * num_lats * num_lons * size_of_float32 * MB_per_bit
# target_chunks_size_MB = int(chnk_sz.split("MB")[0])
# num_chunks = total_size_MB / target_chunks_size_MB
# chnks_per_dim = np.sqrt(num_chunks)
# chnk_lat = int(round(num_lats / chnks_per_dim))
# chnk_lon = int(round(num_lons / chnks_per_dim))
# chnk_time = int(round(tsteps / num_chunks))

for i, f in enumerate(files):
    # open dataset
    if ".grib2" in f:
        if i == 0: # do everything for the grib file in a single iteration of the loop
            path = Path(f)
            fldr_parent = path.parent.absolute()
            year = extract_file_timestep(f).date().year
            if year < 2012:
                origin = "nssl_preciprate_reanalysis_dataset_20010101_to_20111231"
                tsteps = 288 # num of 5 minute timesteps in 24h
            else:
                origin = "iowa_state_mesonet_preciprate_dataset_20150101_to_present"
                tsteps = 720 # num of 2 minute timesteps in 24h

            # create temporary grib combining all the files for a single day
            # delete temporary grib file if it already exists
            try:
                os.remove('{}{}.grib2'.format(fldr_out_tmp_grib, in_date))
            except:
                # file doesn't exist
                pass
            # concatenate grib files
            try:
                output_file = f'{fldr_out_tmp_grib}{in_date}.grib2'
                # Ensure the directory exists
                os.makedirs(os.path.dirname(output_file), exist_ok=True)
                for f in files:
                    with open(output_file, "ab") as myfile, open(f, "rb") as file2:
                        myfile.write(file2.read())
            except:
                sys.exit("Script failed when attempting to concatenat grib files. The number of files on this day was {}. The first file was {} and the last was {}.".format(len(files), files[0], files[-1]))
            # open concatenated grip file and remove unnecessary variables and coordinates
            try:
                ds = xr.open_dataset('{}{}.grib2'.format(fldr_out_tmp_grib, in_date), engine="cfgrib", chunks={"longitude":'auto', "latitude":'auto'},
                                    backend_kwargs={'indexpath': ''})
                try:
                    ds = ds.drop_vars("surface")
                except:
                    pass
                try:
                    del ds.coords['heightAboveSea']
                except:
                    pass
                try:
                    del ds.coords['valid_time']
                except:
                    pass
                try:
                    del ds.coords['step']
                except:
                    pass
            except:
                sys.exit("Failed to create netcdf for {}. Failed loading the consolidated .grib2 data. Data is likely missing on this day.".format(in_date))
            # rename rainfall variable for RainyDay
            ds = ds.rename({"unknown":"rainrate"})
            # change datatypes from float64 to float 32 to save memory
            ds["longitude"] = ds["longitude"].astype(np.float32)
            ds["latitude"] = ds["latitude"].astype(np.float32)
            # define attributes
            ds.attrs['rainrate_units'] = "mm_per_hour"
            ds.attrs['time_units'] = "proceeding"
            ds.attrs['time_zone'] = "UTC"
            ds.attrs['longitude_units'] = "degrees_east"
            ds.attrs['latitude_units'] = "degrees_north"
            ds.attrs['gridcell_feature_represented_by_coordinate'] = "center"
            ds.attrs['source'] = {"origin":origin,
                                "units_information":"https://www.nssl.noaa.gov/projects/mrms/operational/tables.php",
                                "download_page_for_data_from_2001-2011_reanalysis_dataset":"https://edc.occ-data.org/nexrad/mosaic/",
                                "download_page_example_for_data_from_2015_to_present":"https://mtarchive.geol.iastate.edu/2020/09/09/mrms/ncep/PrecipRate/",
                                "note":"""Products at the Iowa State Environmental Mesonet were based on MRMS
            version 10.5-11.5 from 10/2014-9/2020 and version 12-12.1 
            from 10/2020 to present according to Jian Zhang, the first 
            author on the NSSL reanalysis dataset. They are working on 
            a reanalysis dataset from 2012-present which should be ready 
            in 2022 or 2023."""}
            ds.attrs['warnings'] = dic_warnings
        else:
            continue
    else: # this section deals with the quantized netcdfs that actually turned out to be garbage data
        sys.exit("Assuming this is a quantized netcdf, these turned out to be gargage data so not using.")
        # try:
        #     ds = xr.open_dataset(f, chunks={"longitude":'auto', "latitude":'auto'})
        # except:
        #     continue
        # # define attributes
        # ds.attrs['rainrate_units'] = "mm"
        # ds.attrs['time_units'] = "proceeding"
        # ds.attrs['time_zone'] = "UTC"
        # ds.attrs['longitude_units'] = "degrees_west"
        # ds.attrs['latitude_units'] = "degrees_north"
        # ds.attrs['gridcell_feature_represented_by_coordinate'] = "lower_left"
        # ds.attrs['source'] = {"origin":"iowa_state_mesonet_quantized_pngs_to_netcdfs_dataset_20130801_to_20141231",
        #                       "units_information":"https://mesonet.agron.iastate.edu/GIS/rasters.php?rid=3",
        #                       "example_download_link":"https://mesonet.agron.iastate.edu/cgi-bin/request/raster2netcdf.py?dstr=201410250000&prod=mrms_a2m"}
        # # standardize dimensions
        # tstamp = extract_file_timestep(f)
        # ds.coords["time"] = tstamp
        # ds = ds.expand_dims({"time":1})
        # ds = ds.rename({"lon":"longitude", "lat":"latitude", "mrms_a2m":"rainrate"})
        # ds["rainrate"] = ds["rainrate"].astype(np.float32)
        # ds = ds.chunk(chunks={"longitude":'auto', "latitude":'auto', "time":1})   
        # # extract lat and long for first dataset to compare with all the others
        # if i == 0:       
        #     # extract latitude and longitude coordinates
        #     longitude_first_tstep = ds["longitude"].values
        #     latitude_first_tstep = ds["latitude"].values
        # # check to make sure the longitude and latitude values are unchanged
        # else:
        #     long_dif = np.unique(longitude_first_tstep - ds["longitude"].values)
        #     lat_dif = np.unique(latitude_first_tstep - ds["latitude"].values)
        #     if len(long_dif) > 1 or len(lat_dif) > 1:
        #         dic_warnings["grid_is_unchanged_between_timesteps"] = False
        #     if len(ds["longitude"].values) != 7000 or len(ds["latitude"].values) != 3500:
        #         dic_warnings["data_conforms_to_standard_grid_dimensions"] = False
        # ds.attrs['warnings'] = dic_warnings
        # lst_ds.append(ds)

if len(lst_ds) > 0: 
    ds_comb = xr.concat(lst_ds, dim="time") 
else:
    ds_comb = ds

#%% export combined dataset as intermediate output
fl_out_zar = fldr_out_zar_day+"{}.zarr".format(in_date)
print("exporting to zarr....")
bm_time = time.time()
ds_comb = ds_comb.chunk(chunks={"longitude":'auto', "latitude":'auto', "time":'auto'}).to_zarr(fl_out_zar, mode="w")
print(f"successfully exported zarr after time {(time.time() - bm_time)/60:.2f}")
print("Processing dataset for export....")
ds_comb = xr.open_zarr(store=fl_out_zar, chunks={"longitude":'auto', "latitude":'auto', "time":'auto'})

#%% Make final adjustments and quality checks to resulting daily dataset
# sort
ds_comb = ds_comb.sortby(["time", "latitude", "longitude"])

# check timestep regularity:
## infer timestep
v_time = np.sort(ds_comb.time.values)
v_time_diff = np.diff(v_time) / np.timedelta64(1, 'm')
try:
    s_tstep = stats.mode(v_time_diff, axis=0, keepdims=True) # time step scalar
except:
    s_tstep = stats.mode(v_time_diff, axis=0) 
s_tstep = pd.Timedelta(int(s_tstep[0]), 'm')
s_tstep_previous = s_tstep
# check that timestep is standard 2 minutes or 5 minutes; if not, assign it one
ds_comb.warnings["standard_timestep"] = True
if not (s_tstep == pd.Timedelta(2, 'm') or s_tstep == pd.Timedelta(5, 'm')):
    v_time = pd.DatetimeIndex(np.sort(ds_comb.time.values))
    if sum(v_time.minute.values % 2) == 0: # if all the minutes are divisible by 2
        s_tstep = pd.Timedelta(2, 'm')
    elif sum(v_time.minute.values % 5) == 0: # if all the minutes are divisible by 5
        s_tstep = pd.Timedelta(5, 'm')
    else:
        ds_comb.warnings["standard_timestep"] = False
    
ds_comb.attrs['time_step_min']  = str(s_tstep)
ds_comb.attrs['time_step_original_min']  = str(s_tstep_previous)

# check whether the entire 24-hour period is accounted for
v_time_diff = np.append(v_time_diff, s_tstep_previous / np.timedelta64(1, 'm'))
df_counts = pd.DataFrame(np.unique(v_time_diff, return_counts=True)).T
df_counts.columns=["timestep_min", "count"]
# df_counts.astype(int)

s_hours_accounted_for = df_counts.product(axis="columns").sum()/60

if len(df_counts)>1:
    ds_comb.warnings["data_complete_ie_no_missing_timesteps"] = False
if s_hours_accounted_for != 24:
    ds_comb.warnings["data_complete_ie_no_missing_timesteps"] = False

# ensure a proceeding interval
if ds_comb.attrs['time_units'] != "proceeding":
    if ds_comb.attrs['time_units']  == "preceding":
        ds_comb['time'] = ds_comb['time'] - s_tstep
        ds_comb.attrs['time_units']  = "proceeding"
    else:
        lst_problems.append("Time units are unrecognized so conversion to a preceding could not be completed.")
        # sys.exit("Time units are unrecognized so conversion to a preceding could not be completed.")

# # reindex dataset so that all 24 hours are accounted for
if ds_comb.warnings["data_complete_ie_no_missing_timesteps"] == False:
    v_time = pd.DatetimeIndex(np.sort(ds_comb.time.values))
    s_date = (pd.to_datetime(min(v_time))).date()
    try:
        v_time_complete_24h = pd.date_range(start=s_date, end=s_date + pd.Timedelta(1, "day"),\
                                            freq=s_tstep, inclusive='left')
    except:
        v_time_complete_24h = pd.date_range(start=s_date, end=s_date + pd.Timedelta(1, "day"),\
                                    freq=s_tstep, closed='left')
    missing_tsteps = []
    for d in v_time_complete_24h:
        if d not in v_time:
            missing_tsteps.append(str(d))
    #
    ds_comb = ds_comb.reindex(indexers={"time":v_time_complete_24h})
    ds_comb = ds_comb.unify_chunks()
    #
    ds_comb.attrs["missing_timesteps_filled_with_NA_values"] = missing_tsteps
    
# ensure consistent coordinates and dimensions
## rainrate
if ds_comb.attrs['rainrate_units'] != "mm_per_hour":
    if ds_comb.attrs['rainrate_units']  == "mm":
        ds_comb['rainrate'] = ds_comb['rainrate'] * (np.timedelta64(1, "h") / s_tstep) # convert rainrate from mm per timestep to mm per hour by multiplying rainrate by timesteps per hour
        ds_comb.attrs['rainrate_units']  = "mm_per_hour"
    else:
        lst_problems.append("units of measure are unrecognized so conversion to mm per hour could not be completed.")
        # sys.exit("units of measure are unrecognized so conversion to mm per hour could not be completed.")
## longitude
if ds_comb.attrs['longitude_units'] != "degrees_east":
    if ds_comb.attrs['longitude_units']  == "degrees_west":
        ds_comb['longitude'] = 360 +  ds_comb['longitude']
        ds_comb.attrs['longitude_units']  = "degrees_east"
    else:
        lst_problems.append("Longitude units are unrecognized so conversion to degrees east could not be completed.")
        # sys.exit("Longitude units are unrecognized so conversion to degrees east could not be completed.")

if (np.min(ds_comb['longitude'].values) < 0) or (np.max(ds_comb['longitude'].values) > 360):
    lst_problems.append("Longitude values are messed up. This indicates inconsistency in the grid between the netcdfs from quantized pngs or the grib files.")
    # sys.exit("Longitude values are messed up. This indicates inconsistency in the grid between the netcdfs from quantized pngs or the grib files.")

## ensure gridcells represent the center
ds_comb.attrs['center_of_northwest_cell_longitude_deg_east'] = s_nw_corner_lon
ds_comb.attrs['center_of_northwest_cell_latitude_deg_north'] = s_nw_corner_lat
ds_comb.attrs['center_of_southeast_cell_longitude_deg_east'] = s_se_corner_lon
ds_comb.attrs['center_of_southeast_cell_latitude_deg_north'] = s_se_corner_lat
ds_comb.attrs['gricell_center_source'] = "Personal communications with Dr. Jian Zhang, lead on the 2001-2011 MRMS reanalysis dataset (https://doi.org/10.25638/EDC.PRECIP.0001)" 

def check_for_center(ds_comb, s_nw_corner_lat=s_nw_corner_lat, s_nw_corner_lon=s_nw_corner_lon,
                    s_se_corner_lat=s_se_corner_lat, s_se_corner_lon=s_se_corner_lon):
    # assumes all latitude and longitude values have the same coordinate system
    v_lats = ds_comb.latitude.values
    v_lons = ds_comb.longitude.values
    s_ds_nw_corner_lat = np.round((max(v_lats)), 3)
    s_ds_nw_corner_lon = np.round(min(v_lons), 3)
    s_ds_se_corner_lat = np.round(min(v_lats), 3)
    s_ds_se_corner_lon = np.round(max(v_lons), 3)
    s_ds_minus_center_nw_corner_lat = np.round(s_ds_nw_corner_lat - s_nw_corner_lat, 3)
    s_ds_minus_center_nw_corner_lon = np.round(s_ds_nw_corner_lon - s_nw_corner_lon, 3)
    s_ds_minus_center_se_corner_lat = np.round(s_ds_se_corner_lat - s_se_corner_lat, 3)
    s_ds_minus_center_se_corner_lon = np.round(s_ds_se_corner_lon - s_se_corner_lon, 3)
    return s_ds_minus_center_nw_corner_lat, s_ds_minus_center_nw_corner_lon, s_ds_minus_center_se_corner_lat, s_ds_minus_center_se_corner_lon

if ds_comb.attrs['gridcell_feature_represented_by_coordinate'] != "center":
    if ds_comb.attrs['gridcell_feature_represented_by_coordinate'] == "lower_left":
        # extract current latitude and longitude values
        v_lats = ds_comb.latitude.values
        v_lons = ds_comb.longitude.values
        # shift coordinates east 1/2 gridcell
        v_lats_diff = abs(np.diff(v_lats))
        v_lats_diff = np.append(v_lats_diff, v_lats_diff[-1])
        v_lats_correction = v_lats_diff/2
        v_lats = v_lats + v_lats_correction
        #
        # shift coordinates north 1/2 gridcell
        v_lon_diff = abs(np.diff(v_lons))
        v_lon_diff = np.append(v_lon_diff, v_lon_diff[-1])
        v_lon_correction = v_lon_diff/2
        v_lons = v_lons + v_lon_correction
        #
        # update the latitude and longitude coordinates
        ds_comb["latitude"] = v_lats 
        ds_comb["longitude"] = v_lons
        #
        ds_comb.attrs['gridcell_feature_represented_by_coordinate'] = "center"
    else:
        lst_problems.append("The gricell feature represented by the coordinates is neither 'center' nor 'lower_left'")

### if the coordinates are classified as representing the center but the sum of the differences is not 0:
s_ds_minus_center_nw_corner_lat, s_ds_minus_center_nw_corner_lon, s_ds_minus_center_se_corner_lat, s_ds_minus_center_se_corner_lon = check_for_center(ds_comb)
if (abs(s_ds_minus_center_nw_corner_lat) + abs(s_ds_minus_center_nw_corner_lon) + \
        abs(s_ds_minus_center_se_corner_lat) + abs(s_ds_minus_center_se_corner_lon)) != 0:
    # if these sum to anything other than 0
    s_problem = """The gridcell feature represented by coordinates is set to 'center' 
    but there is a discrepancy between either the northwestern-most or southeastern-most coordinates.
    Northwest lat and lon differences: ({}, {})
    Southeast lat and lon differences: ({}, {})
    """.format(s_ds_minus_center_nw_corner_lat, s_ds_minus_center_nw_corner_lon, 
        s_ds_minus_center_se_corner_lat, s_ds_minus_center_se_corner_lon)
    lst_problems.append(s_problem)

## Check to see if the lats and lons conform to the most recent grid and 
## assign it the most recent grid if not
v_lats = ds_comb.latitude.values
v_lons = ds_comb.longitude.values
chk_lats = sum(v_lats != v_lats_most_recent.values)
chk_lons = sum(v_lons != v_lons_most_recent.values)
s_mean_lat_dif = np.mean(v_lats_most_recent.values - v_lats)
s_mean_lon_dif = np.mean(v_lons_most_recent.values - v_lons)
ds_comb.attrs["coordinates_adjusted_to_conform_to_modern_grid"] = "False"
if chk_lats>0 or chk_lons>0:
    ds_comb.attrs["coordinates_adjusted_to_conform_to_modern_grid"] = "True"
    ds_comb.attrs["previous_latitude_values_on_grid_centers"] = v_lats.tolist()
    ds_comb.attrs["previous_longitude_values_on_grid_centers"] = v_lons.tolist()
    ds_comb.attrs["mean_latitude_shift_in_translation_to_modern_grid (most recent minus previous, so positive means the most recent grid was more north)"] = float(s_mean_lat_dif)
    ds_comb.attrs["mean_longitude_shift_in_translation_to_modern_grid (most recent minus previous, so positive means the most recent grid was more east)"] = float(s_mean_lon_dif)
    ds_comb["latitude"] = v_lats_most_recent
    ds_comb["longitude"] = v_lons_most_recent

## Convert to upper left coordinates
### extract current latitude and longitude values
v_lats = ds_comb.latitude.values
v_lons = ds_comb.longitude.values
### shift coordinates north 1/2 gridcell
v_lats_diff = abs(np.diff(v_lats))
v_lats_diff = np.append(v_lats_diff, v_lats_diff[-1])
v_lats_correction = v_lats_diff/2
v_lats = v_lats + v_lats_correction

### shift coordinates west 1/2 gridcell
v_lon_diff = abs(np.diff(v_lons))
v_lon_diff = np.append(v_lon_diff, v_lon_diff[-1])
v_lon_correction = v_lon_diff/2
v_lons = v_lons - v_lon_correction

### update the latitude and longitude coordinates
ds_comb["latitude"] = v_lats 
ds_comb["longitude"] = v_lons

ds_comb.attrs['gridcell_feature_represented_by_coordinate'] = "upper_left"
ds_comb.attrs['grid_spacing'] = np.mean(v_lon_diff).round(3)

## convert any negative values to nan
ds_comb = ds_comb.where(ds_comb["rainrate"]>=0) # converts "falses" to nan, where false means less than 0

## convert any values above 9000 to nan
ds_comb = ds_comb.where(ds_comb["rainrate"]<9000) # converts "falses" to nan, where false means greater than 9000

## Append attributes created in this script
ds_comb.attrs["problems"] = "none detected" if len(lst_problems)==0 else str(lst_problems)
ds_comb.attrs["source"] = str(ds_comb.attrs["source"])
# create an attribute for each warning
for key in ds_comb.attrs["warnings"]:
    ds_comb.attrs[key] = ds_comb.attrs["warnings"][key]
del ds_comb.attrs["warnings"]
# ds_comb.attrs["warnings"] = str(ds_comb.attrs["warnings"])


## update attributes
ds_comb.latitude.attrs["units"] = "degrees_north"
ds_comb.longitude.attrs["units"] = "degrees_east"

### clear out rainrate attributes
current_attributes = ds_comb.rainrate.attrs.copy()
for at in current_attributes:
    try:
        del ds_comb.rainrate.attrs[at]
    except:
        continue

ds_comb.rainrate.attrs["long_name"] = "Precipitation Rate"
ds_comb.rainrate.attrs["short_name"] = "PrecipRate"
ds_comb.rainrate.attrs["units"] = "mm/hr"
ds_comb.rainrate.attrs["description"] = "Radar Precipitation Rate"
ds_comb.rainrate.attrs["missing_value"] = np.nan
ds_comb.rainrate.attrs["_FillValue"] = np.nan
ds_comb.rainrate.attrs["Grib2_Parameter_Name"] = "PrecipRate"
ds_comb.rainrate.attrs["Originating_or_generating_Center"] = "US NOAA Office of Oceanic and Atmospheric Research"

# print("Loaded and formatted dataset from raw data: {}".format(time.time() - bm_time))
#%% export to netcdf
#%% first write a zarr file (only necessary for the netcdf files)
print(f"Prepared dataset for exporting. Time elapsed: {(time.time() - start_time)/60:.2f}")
# remove unnecesary coordinates and attributes
ds_comb = remove_vars(ds_comb)
bm_time = time.time()
ds_comb.to_netcdf(fl_out_nc, encoding= {"rainrate":{"zlib":True}})
print(f"successfully exported netcdf after time {(time.time() - bm_time)/60:.2f}")
bm_time = time.time()
shutil.rmtree(fl_out_zar)
if ".grib2" in f:
    # delete temporary concatenated grib file
    os.remove('{}{}.grib2'.format(fldr_out_tmp_grib, in_date))
#%% final benchmark
elapsed = time.time() - start_time
print("Succeeded in creating netcdf for {}. Total script runtime: {:.2f}".format(in_date, elapsed/60))
# print("Total script run time: {} seconds.".format(elapsed, in_date))