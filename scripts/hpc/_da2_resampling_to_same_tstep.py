#%% Import libraries
import time
import shutil
import xarray as xr
import geopandas as gp
import pandas as pd
import zarr
import sys
from pathlib import Path
import gc
import dask
dask.config.set(**{'array.slicing.split_large_chunks': False}) # to silence warnings of loading large slice into memory
# dask.config.set(scheduler='synchronous') # this forces single threaded computations
dask.config.set(scheduler='threads')
from __utils import *
import ast
from glob import glob
from rasterio.enums import Resampling
import rioxarray
import numpy as np

import warnings
warnings.filterwarnings("ignore")

start_time = time.time()

overwrite_existing_outputs = True
bias_correction_reference = "aorc"
target_tstep_min = 2

tsteps_per_day = int(24 * 60 / target_tstep_min)
tsteps_per_hr = 60 / target_tstep_min
final_chunking_dict = dict(time = tsteps_per_day/2, latitude = 1, longitude = 1)
# WITH exporting penultimate zarr
# final_chunking_dict = dict(time = tsteps_per_day/2, latitude = 5, longitude = 5)
# final_chunking_dict = dict(time = tsteps_per_day/2, latitude = 3, longitude = 3)
# final_chunking_dict = dict(time = tsteps_per_day/2, latitude = 1, longitude = 1)
# final_chunking_dict = dict(time = tsteps_per_day, latitude = 4, longitude = 4) | killed
# final_chunking_dict = dict(time = tsteps_per_day, latitude = 3, longitude = 3) | killed
# final_chunking_dict = dict(time = tsteps_per_day, latitude = 1, longitude = 1) | killed
# final_chunking_dict = dict(time = tsteps_per_day, latitude = 500, longitude = 500) | killed
# final_chunking_dict = dict(time = tsteps_per_day, latitude = 50, longitude = 50) | killed
# final_chunking_dict = dict(time = tsteps_per_day, latitude = 5, longitude = 5) | killed
# WITHOUT exporting penultimate zarr
# final_chunking_dict = dict(time = tsteps_per_day, latitude = 1, longitude = 1) | killed
# final_chunking_dict = dict(time = tsteps_per_hr*2, latitude = 500, longitude = 500) | killed
# final_chunking_dict = dict(time = tsteps_per_day/2, latitude = 10, longitude = 10) | killed
# final_chunking_dict = dict(time = tsteps_per_hr, latitude = 100, longitude = 100) | script finished. Elapsed time (min): 45.73
# final_chunking_dict = dict(time = tsteps_per_hr, latitude = 500, longitude = 500) | script finished. Elapsed time (min): 49.03
# final_chunking_dict = dict(time = tsteps_per_day/2, latitude = 50, longitude = 50) | killed
# final_chunking_dict = dict(time = tsteps_per_day, latitude = 4, longitude = 4) | killed
# final_chunking_dict = dict(time = tsteps_per_day, latitude = 5, longitude = 5) | killed
# final_chunking_dict = dict(time = 1, latitude = 3500, longitude = 3500) | script finished. Elapsed time (min): 50.91
# final_chunking_dict = dict(time = tsteps_per_day, latitude = 50, longitude = 50) | killed
# final_chunking_dict = dict(time = tsteps_per_day, latitude = 500, longitude = 500) | killed

final_output_type = "zarr" # must be "nc" or "zarr"

# chnk_sz = "100MB"

performance = {}
#%% work

# in_date = "20210721" # "20210719" corresponds to slurm task array 200
# fldr_zarr_fullres_daily = "/project/quinnlab/dcl3nd/norfolk/highres-radar-rainfall-processing/data/mrms_zarr_preciprate_fullres_dailyfiles/"
# fldr_zarr_fullres_daily_constant_tstep = "/project/quinnlab/dcl3nd/norfolk/highres-radar-rainfall-processing/data/mrms_zarr_preciprate_fullres_dailyfiles_constant_tstep/"
# fldr_scratch_zarr = "/scratch/dcl3nd/highres-radar-rainfall-processing/_scratch/zarrs/"
# fldr_scratch_csv = "/scratch/dcl3nd/highres-radar-rainfall-processing/_scratch/csv/"
# fldr_bias_crxn_ref = "/project/quinnlab/dcl3nd/norfolk/highres-radar-rainfall-processing/data/raw_data/raw_data/aorc/"
# # fldr_bias_crxn_ref = "/project/quinnlab/dcl3nd/norfolk/stormy/data/climate/StageIV_rainfall/"
# # f_shp_sst_transom = "/project/quinnlab/dcl3nd/stormy/stochastic_storm_transposition/norfolk/transposition_domain/norfolk_trans_dom_4326.shp"
# f_shp_sst_transom = None
#%% end work

# folders (with proceeding fwd slash)
in_date = str(sys.argv[1]) # YYYYMMDD
fldr_zarr_fullres_daily = str(sys.argv[2]) # ${assar_dirs[out_fullres_dailyfiles]} # "/project/quinnlab/dcl3nd/norfolk/highres-radar-rainfall-processing/data/mrms_zarr_preciprate_fullres_dailyfiles/"
fldr_zarr_fullres_daily_constant_tstep = str(sys.argv[3]) # ${assar_dirs[out_fullres_dailyfiles_consolidated]} # "/project/quinnlab/dcl3nd/highres-radar-rainfall-processing/mrms_zarr_preciprate_fullres_dailyfiles_constant_tstep/"
fldr_scratch_zarr = str(sys.argv[4]) # ${assar_dirs[scratch_zarrs]} # "/scratch/dcl3nd/highres-radar-rainfall-processing/_scratch/zarrs/"
fldr_scratch_csv = str(sys.argv[5]) # ${assar_dirs[scratch_zarrs]} # "/project/quinnlab/dcl3nd/norfolk/highres-radar-rainfall-processing/data/_scratch/csv/"
fldr_bias_crxn_ref = str(sys.argv[6])

try:
    f_shp_sst_transom = str(sys.argv[7]) # ${assar_dirs[shp_transposition_domain]} # "/project/quinnlab/dcl3nd/norfolk/stormy/stochastic_storm_transposition/norfolk/transposition_domain/norfolk_trans_dom_4326.shp"
except:
    f_shp_sst_transom = None
    print("Bias correcting entire MRMS dataset since no domain shapefile has been supplied........")
performance["date"] = in_date
print(f"Processing date {in_date}")
# f_out_export_perf = fldr_scratch_zarr + "_export_stats_{}.csv".format(in_date)
#%%  
fl_in_zarr = fldr_zarr_fullres_daily +"{}.zarr".format(in_date)
fl_out_zarr = fldr_zarr_fullres_daily_constant_tstep +"{}.zarr".format(in_date)
fl_out_csv = fldr_scratch_csv +"da2_resampling_{}.csv".format(in_date)
fl_out_csv_qaqc = fldr_scratch_csv +"qaqc_of_daily_fullres_data_{}.csv".format(in_date)

# if (not overwrite_existing_outputs) and Path(fl_out_zarr).exists():
#     print(f"File already exists and overwrite_existing_outputs is set to {overwrite_existing_outputs}. Not reprocessing {fl_out_zarr}")
#     sys.exit(0)

if not Path(fl_in_zarr).exists():
    print(f"Raw zarr file not found. Skipping {fl_in_zarr}")
    sys.exit(0)

# find related stageIV rainfall file

if bias_correction_reference == "stageiv":
    f_nc_stageiv = glob(fldr_bias_crxn_ref + "{}/*{}*".format(in_date[0:4],in_date))
    if len(f_nc_stageiv)>0:
        f_nc_stageiv = f_nc_stageiv[0]
        # data_available_for_bias_correction = True
        performance["data_available_for_bias_correction"] = True
    else:
        # data_available_for_bias_correction = False
        performance["data_available_for_bias_correction"] = False
        print("No stage iv data for this date. No bias correction being performed....")

if bias_correction_reference == "aorc":
    year=in_date[0:4]
    lst_f_aorc_yr = glob(f"{fldr_bias_crxn_ref}/data/{year}*")
    if len(lst_f_aorc_yr) == 0:
        performance["data_available_for_bias_correction"] = False
        # data_available_for_bias_correction = False
    else:
        performance["data_available_for_bias_correction"] = True
        # data_available_for_bias_correction = True
        f_aorc_yr = lst_f_aorc_yr[0]
        da_aorc_rainfall = xr.open_dataset(f_aorc_yr, engine = "zarr", chunks = 'auto')["APCP_surface"]
        da_aorc_rainfall = da_aorc_rainfall.sel(time=in_date)
# if toggled to not reprocess existing outputs, only continue if they don't already exist
if not overwrite_existing_outputs:
    try:
        s_perf = pd.read_csv(fl_out_csv).iloc[0]
        if final_output_type == "nc":
            check = (s_perf["to_nc_errors"] == "None") and (s_perf["problem_exporting_nc"] == False)
        elif final_output_type == "zarr":
            check = (s_perf["to_zarr_errors"] == "None") and (s_perf["problem_exporting_zarr"] == False)
        if check:
            print(f"{fl_out_zarr} successfully generated with no errors and overwrite_existing_outputs is set to {overwrite_existing_outputs}. Not re-running.")
            sys.exit(0)
    except:
        pass

if f_shp_sst_transom is not None:
    gdf_transdomain = gp.read_file(f_shp_sst_transom)
else: 
    gdf_transdomain = None

#%% functions
def define_zarr_compression(ds, clevel=5):
    encoding = {}
    for da_name in ds.data_vars:
        encoding[da_name] = {"compressor": zarr.Blosc(cname="zstd", clevel=clevel, shuffle=zarr.Blosc.SHUFFLE)}
    return encoding

def clip_ds_to_transposition_domain(ds, gdf_transdomain, buffer = 0.15):
    # clips a rectangle that overs the transposition domain
    transdom_bounds = gdf_transdomain.bounds
    ds = ds.where((ds.latitude >= float(transdom_bounds.miny.iloc[0]-.05-buffer)) & (ds.latitude <= float(transdom_bounds.maxy.iloc[0]+.05+buffer)) & (ds.longitude >= float(transdom_bounds.minx.iloc[0]+360-.05-buffer)) & (ds.longitude <= float(transdom_bounds.maxx.iloc[0]+360+.05+buffer)), drop = True)
    return ds

def compute_total_rainfall_over_domain(ds):
    # find daily total rainfall over entire domain in mm
    # assumes hourly rainfall intensities
    tot_rain = ds.rainrate.mean(dim = ["time", "latitude", "longitude"])*24
    return tot_rain.values

def estimate_chunk_memory(ds, input_chunk_sizes=None):
    """
    Estimate memory requirements for each chunk of an Xarray dataset.

    Parameters:
    -----------
    ds : xr.Dataset or xr.DataArray
        The Xarray dataset or data array.
    chunk_sizes : dict, optional
        A dictionary specifying the chunk sizes for each dimension (e.g., {'time': 100, 'latitude': 500}).
        If None, it uses the current chunk sizes in the dataset.

    Returns:
    --------
    float
        The estimated memory usage per chunk in megabytes (MB).
    """
    import numpy as np
    #
    # Use existing chunk sizes if none are provided
    if input_chunk_sizes is None:
        chunk_sizes = {dim: ds.chunks.get(dim, (len(ds[dim]),))[0] for dim in ds.dims}
    else:
        chunk_sizes = input_chunk_sizes.copy()
        # assume full chunking on all other datasets
        keys_to_skip = chunk_sizes.keys()
        for dim in ds.dims:
            if dim not in keys_to_skip:
                chunk_sizes[dim] = len(ds[dim])
    #
    # Estimate the total number of elements in one chunk
    total_elements = np.prod(list(chunk_sizes.values()))
    #
    # Get the dtype of the dataset (assuming homogeneous dtype across variables)
    dtype = ds.dtype if isinstance(ds, xr.DataArray) else next(iter(ds.data_vars.values())).dtype
    #
    # Calculate the size of each element in bytes
    element_size = np.dtype(dtype).itemsize
    # Estimate the total memory for one chunk (in bytes)
    total_bytes = total_elements * element_size
    # Convert to megabytes (MB) for convenience
    total_mb = total_bytes / (1024 ** 2)
    # print(f"Estimated memory per chunk: {total_mb:.2f} MB")
    return total_mb, chunk_sizes

def bias_correct_and_fill_mrms(ds_mrms, ds_ref, lst_tmp_files_to_delete, 
                               crxn_upper_bound = crxn_upper_bound, crxn_lower_bound = crxn_lower_bound,
                                 verbose = False):
    mrms_time_encoding = {k: ds_mrms.time.encoding[k] for k in ['units', 'calendar', 'dtype'] if k in ds_mrms.time.encoding}
    # convert stage iv to proceeding time interval
    # ds_ref_og = xr.open_dataset(tmp_raw_ref_zarr, chunks = dic_auto_chunk, engine = "zarr")
    # ds_ref = ds_ref_og.copy()
    # create hourly version of mrms data
    try:
        ds_mrms_hourly = ds_mrms.resample(time = "H").mean() # convert to hourly timestep
    except Exception as e:
        print("there was an errro when computing this line:")
        print("ds_mrms_hourly = ds_mrms.resample(time = 'H').mean()")
        print(e)
        print(f"MRMS dims: {ds_mrms.dims}")
        print(f"MRMS coords: {ds_mrms.coords}")
        sys.exit()
    # spatially resample MRMS data to match stage IV resolution
    total_mb_ref, dic_chunks_ref = estimate_chunk_memory(ds_ref)
    total_mb_mrms, dic_chunks_mrms = estimate_chunk_memory(ds_mrms)
    xds_mrms_hourly_to_ref = spatial_resampling(ds_mrms_hourly.chunk(dic_chunks_ref), ds_ref, "latitude", "longitude")
    #
    ## write the bias correction dataset to a temporary file
    # tmp_zarr = f"{fldr_scratch_zarr}{in_date}_xds_mrms_hourly_to_ref.zarr"
    # lst_tmp_files_to_delete.append(tmp_zarr)
    # # gc.collect()
    # bm_time = time.time()
    # xds_mrms_hourly_to_ref.chunk(dict(time = "auto", latitude = "auto", longitude = "auto")).to_zarr(tmp_zarr, mode = "w", encoding = define_zarr_compression(xds_mrms_hourly_to_ref))
    # xds_mrms_hourly_to_ref = xr.open_zarr(store=tmp_zarr).chunk(dict(time = "auto", latitude = "auto", longitude = "auto"))
    # # gc.collect()
    # # print("exported xds_mrms_hourly_to_ref to zarr")
    # print(f"Time to export xds_mrms_hourly_to_ref to zarr (min): {((time.time() - bm_time)/60):.2f} | total script runtime (min): {((time.time() - start_time)/60):.2f}")
    #
    #
    # compute correction factor
    xds_mrms_hourly_correction_factor_refres = ds_ref/xds_mrms_hourly_to_ref
    # if mrms is 0, assign a 1 correction factor (the case where stage iv is non zero is taken care of below) (otherwise it will be infinity)
    xds_mrms_hourly_correction_factor_refres = xr.where((xds_mrms_hourly_to_ref <= 0),
                                                        x = 1, y = xds_mrms_hourly_correction_factor_refres)
    # to mitigate outliers, enforce correction factor bounds
    ## where stage iv is 0 or negative, assign a value of 1 (so no correction, assuming stage iv is missing in these locations)
    ### this should only affect zeros since I'm covnerting negative values to 0
    xds_mrms_hourly_correction_factor_refres = xr.where((ds_ref <= 0),
                                                        x = 1, y = xds_mrms_hourly_correction_factor_refres)
    ## where upper bound is exceeded, assign upper bound
    xds_mrms_hourly_correction_factor_refres = xr.where((xds_mrms_hourly_correction_factor_refres > crxn_upper_bound),
                                                        x = crxn_upper_bound, y = xds_mrms_hourly_correction_factor_refres)
    ## where lower bound is exceeded, assign lower bound
    xds_mrms_hourly_correction_factor_refres = xr.where((xds_mrms_hourly_correction_factor_refres < crxn_lower_bound),
                                                        x = crxn_lower_bound, y = xds_mrms_hourly_correction_factor_refres)
    # upsample correction factor in space to MRMS resolution
    xds_mrms_hourly_correction_factor_fulres= spatial_resampling(xds_mrms_hourly_correction_factor_refres, ds_mrms_hourly, "latitude", "longitude")
    #
    ## write the bias correction dataset to a temporary file
    # tmp_zarr = f"{fldr_scratch_zarr}{in_date}_xds_mrms_hourly_correction_factor_fulres.zarr"
    # lst_tmp_files_to_delete.append(tmp_zarr)
    # bm_time = time.time()
    # print("exporting xds_mrms_hourly_correction_factor_fulres with chunk size and chunks: {}, {}".format(total_mb_mrms, dic_chunks_mrms))
    # # gc.collect()
    # xds_mrms_hourly_correction_factor_fulres.chunk(dict(dic_chunks_mrms)).to_zarr(tmp_zarr, mode = "w", encoding = define_zarr_compression(xds_mrms_hourly_correction_factor_fulres))
    # xds_mrms_hourly_correction_factor_fulres = xr.open_zarr(store=tmp_zarr).chunk(dict(dic_chunks_mrms))
    # # gc.collect()
    # print("exported xds_mrms_hourly_correction_factor_fulres to zarr")
    # print(f"Time to export (min): {((time.time() - bm_time)/60):.2f} | total script runtime (min): {((time.time() - start_time)/60):.2f}")
    #
    #
    # if ref is non zero and mrms is zero, fill with stage iv precip intensities
    ### create stage iv data with same res as mrms
    xds_ref_to_mrms= spatial_resampling(ds_ref, ds_mrms_hourly, "latitude", "longitude")
    #
    ## write temporary file
    # bm_time = time.time()
    # tmp_zarr = f"{fldr_scratch_zarr}{in_date}_xds_ref_to_mrms.zarr"
    # lst_tmp_files_to_delete.append(tmp_zarr)
    # # gc.collect()
    # xds_ref_to_mrms.chunk(dict(time = "auto", latitude = "auto", longitude = "auto")).to_zarr(tmp_zarr, mode = "w",
    #                                                                                                                     encoding = define_zarr_compression(xds_ref_to_mrms))
    # xds_ref_to_mrms = xr.open_zarr(store=tmp_zarr).chunk(dict(time = "auto", latitude = "auto", longitude = "auto"))
    # # gc.collect()
    # # print("exported xds_ref_to_mrms to zarr")
    # print(f"Time to export xds_ref_to_mrms to zarr (min): {((time.time() - bm_time)/60):.2f} | total script runtime (min): {((time.time() - start_time)/60):.2f}")
    #
    #
    ### create dataset with ref rainfall intensities at the indieces where condition is true and 0 everywhere else
    xds_ref_where_mrms_is_0_and_ref_is_not_hourly = xr.where((xds_ref_to_mrms > 0) & (ds_mrms_hourly == 0),
                                                            x = xds_ref_to_mrms, y = 0)
    #### upsample to the full MRMS resolution
    ##### forward fill missing values since it is a proceeding dataset
    xds_ref_where_mrms_is_0_and_ref_is_not = xds_ref_where_mrms_is_0_and_ref_is_not_hourly.reindex(dict(time = ds_mrms.time)).ffill(dim="time").chunk(dict(time = "auto", latitude = "auto", longitude = "auto"))
    #
    # print("REMOVED exporting xds_ref_where_mrms_is_0_and_ref_is_not (this is around a common kill point)")
    # bm_time = time.time()
    # encoding = define_zarr_compression(xds_ref_where_mrms_is_0_and_ref_is_not)
    # encoding['time'] = mrms_time_encoding
    # tmp_zarr= f"{fldr_scratch_zarr}{in_date}_xds_ref_where_mrms_is_0_and_ref_is_not.zarr"
    # lst_tmp_files_to_delete.append(tmp_zarr)
    # # gc.collect()
    # xds_ref_where_mrms_is_0_and_ref_is_not.chunk(dict(time = "auto", latitude = "auto", longitude = "auto")).to_zarr(tmp_zarr, mode = "w", encoding = encoding)
    # xds_ref_where_mrms_is_0_and_ref_is_not = xr.open_zarr(store=tmp_zarr).chunk(dict(time = "auto", latitude = "auto", longitude = "auto"))
    # # gc.collect()
    # print("exported xds_ref_where_mrms_is_0_and_ref_is_not to zarr")
    # print(f"Time to export (min): {((time.time() - bm_time)/60):.2f} | total script runtime (min): {((time.time() - start_time)/60):.2f}")
    #
    ### upsample bias correction to full res data
    xds_correction_to_mrms = xds_mrms_hourly_correction_factor_fulres.chunk(dict(dic_chunks_mrms)).reindex(dict(time = ds_mrms.time)).ffill(dim="time")
    #
    # write the bias correction dataset to a temporary file
    tmp_bias_correction_factor = f"{fldr_scratch_zarr}{in_date}_bias_crxn_factor.zarr"
    lst_tmp_files_to_delete.append(tmp_bias_correction_factor)
    # gc.collect()
    # time_before_export = pd.Series(xds_correction_to_mrms.time.values)
    # print("(NECESSARY INTERMEDIATE OUTPUT) exporting xds_correction_to_mrms with chunk size and chunks: {}, {}".format(total_mb_mrms, dic_chunks_mrms))
    bm_time = time.time()
    encoding = define_zarr_compression(xds_correction_to_mrms)
    encoding['time'] = mrms_time_encoding
    # print(f"assigning time encoding to xds_correction_to_mrms before export: {encoding['time']}")
    xds_correction_to_mrms.chunk(dict(dic_chunks_mrms)).to_zarr(tmp_bias_correction_factor, mode = "w", encoding = encoding)
    xds_correction_to_mrms = xr.open_zarr(store=tmp_bias_correction_factor).chunk(dict(dic_chunks_mrms))
    # time_after_export = pd.Series(xds_correction_to_mrms.time.values)
    # print("(necessary) exported bias correction factor to zarr")
    print(f"(necessary) Time to export bias correction factor to zarr (min): {((time.time() - bm_time)/60):.2f} | total script runtime (min): {((time.time() - start_time)/60):.2f}")
    gc.collect()
    #
    ### apply correction factor
    xds_mrms_biascorrected = (ds_mrms * xds_correction_to_mrms)
    #
    # bm_time = time.time()
    # tmp_bias_crctd = f"{fldr_scratch_zarr}{in_date}_bias_crctd.zarr"
    # lst_tmp_files_to_delete.append(tmp_bias_crctd)
    # # gc.collect()
    # xds_mrms_biascorrected.chunk(dict(time = "auto", latitude = "auto", longitude = "auto")).to_zarr(tmp_bias_crctd, mode = "w", encoding = define_zarr_compression(xds_mrms_biascorrected))
    # xds_mrms_biascorrected = xr.open_zarr(store=tmp_bias_crctd).chunk(dict(time = "auto", latitude = "auto", longitude = "auto"))
    # # gc.collect()
    # print("exported xds_mrms_biascorrected to zarr")
    # print(f"Time to export (min): {((time.time() - bm_time)/60):.2f} | total script runtime (min): {((time.time() - start_time)/60):.2f}")
    #
    ### fill in with stageIV data where mrms data is missing
    xds_mrms_biascorrected_filled = xds_mrms_biascorrected + xds_ref_where_mrms_is_0_and_ref_is_not
    #
    bm_time = time.time()
    tmp_bias_crctd_fld = f"{fldr_scratch_zarr}{in_date}_bias_crctd_fld.zarr"
    lst_tmp_files_to_delete.append(tmp_bias_crctd_fld)
    # gc.collect()
    xds_mrms_biascorrected_filled.chunk(dict(time = "auto", latitude = "auto", longitude = "auto")).to_zarr(tmp_bias_crctd_fld, mode = "w", encoding = define_zarr_compression(xds_mrms_biascorrected_filled))
    xds_mrms_biascorrected_filled = xr.open_zarr(store=tmp_bias_crctd_fld).chunk(dict(time = "auto", latitude = "auto", longitude = "auto"))
    # gc.collect()
    # print("(necessary) exported xds_mrms_biascorrected_filled dataset to zarr")
    print(f"(necessary) Time to export xds_mrms_biascorrected_filled dataset to zarr (min): {((time.time() - bm_time)/60):.2f} | total script runtime (min): {((time.time() - start_time)/60):.2f}")
    #
    ### keep original mrms data
    # xds_mrms_biascorrected_filled = xds_mrms_biascorrected_filled.assign(rainrate_uncorrected = ds_mrms.rainrate)
    ### include bias correction ds
    # xds_mrms_biascorrected_filled = xds_mrms_biascorrected_filled.assign(mrms_bias_correction_multiplier = xds_correction_to_mrms.rainrate)
    ### include ref fill values
    # xds_mrms_biascorrected_filled = xds_mrms_biascorrected_filled.assign(ref_fillvals_where_mrms_is_0_and_ref_is_not = xds_ref_where_mrms_is_0_and_ref_is_not.rainrate)
    # check
    ## compare stage iv with bias corrected total rainfall for the whole day 
    ### computing daily totals by finding the average daily intensity in mm per hour with the .mean("time") function and then multiplying by 24 hours
    #### comparing using mrms resolution
    tot_rain_mrms_corrected = compute_total_rainfall_over_domain(xds_mrms_biascorrected_filled.chunk(dict(time = "auto", latitude = "auto", longitude = "auto")))
    tot_rain_ref = compute_total_rainfall_over_domain(ds_ref)
    tot_rain_mrms_uncorrected = compute_total_rainfall_over_domain(ds_mrms_hourly)
    ## old method
    # tot_rain_mrms_corrected = (xds_mrms_biascorrected_filled.mean("time")*24).rainrate.sum().values
    # tot_rain_ref = (xds_ref_to_mrms.mean("time")*24).rainrate.sum().values
    # tot_rain_mrms_uncorrected = (ds_mrms.mean("time")*24).rainrate.sum().values
    # print("Fraction of domain-wide rainfall totals:")
    # print("Bias corrected MRMS data over stage iv data: {}".format(tot_rain_mrms_corrected/tot_rain_ref))
    # print("UNbias corrected MRMS data over stage iv data: {}".format(tot_rain_mrms_uncorrected/tot_rain_ref))
    #### comparing using ref resolution
    ### resample bias corrected mrms data to stage iv resolution to compare
    # xds_mrms_biascorrected_to_ref= spatial_resampling(xds_mrms_biascorrected_filled, ds_ref, "latitude", "longitude")
    # tot_rain_mrms_corrected = (xds_mrms_biascorrected_to_ref.mean("time")*24).rainrate.sum().values
    # tot_rain_mrms_uncorrected = (xds_mrms_biascorrected_to_ref.mean("time")*24).rainrate_uncorrected.sum().values
    # tot_rain_ref = (ds_ref.mean("time")*24).rainrate.sum().values
    # print("Fraction of domain-wide rainfall totals:")
    # print("Bias corrected MRMS data over stage iv data: {}".format(tot_rain_mrms_corrected/tot_rain_ref))
    # print("UNbias corrected MRMS data over stage iv data: {}".format(tot_rain_mrms_uncorrected/tot_rain_ref))
    if tot_rain_ref > 0:
        domainwide_totals_CORRECTED_mrms_over_ref = tot_rain_mrms_corrected/tot_rain_ref
        domainwide_totals_uncorrected_mrms_over_ref = tot_rain_mrms_uncorrected/tot_rain_ref
    else:
        domainwide_totals_CORRECTED_mrms_over_ref = domainwide_totals_uncorrected_mrms_over_ref = np.nan
    xds_mrms_biascorrected_filled.attrs["domainwide_totals_CORRECTED_mrms_mm"] = tot_rain_mrms_corrected
    xds_mrms_biascorrected_filled.attrs["domainwide_totals_ref_mm"] = tot_rain_ref
    xds_mrms_biascorrected_filled.attrs["domainwide_totals_uncorrected_mrms_mm"] = tot_rain_mrms_uncorrected
    xds_mrms_biascorrected_filled.attrs["domainwide_totals_CORRECTED_mrms_over_ref"] = domainwide_totals_CORRECTED_mrms_over_ref
    xds_mrms_biascorrected_filled.attrs["domainwide_totals_uncorrected_mrms_over_ref"] = domainwide_totals_uncorrected_mrms_over_ref
    xds_mrms_biascorrected_filled.attrs["correction_factor_upperbound"] = crxn_upper_bound
    xds_mrms_biascorrected_filled.attrs["correction_factor_lowerbound"] = crxn_lower_bound
    return xds_mrms_biascorrected_filled, xds_mrms_hourly_to_ref, ds_ref,\
            xds_correction_to_mrms, xds_ref_where_mrms_is_0_and_ref_is_not, lst_tmp_files_to_delete

def process_bias_corrected_dataset(ds_mrms_biascorrected_filled, ds_mrms, ds_ref_proceeding, ds_correction_to_mrms,
                                    ds_ref_where_mrms_is_0_and_ref_is_not, lst_tmp_files_to_delete,
                                      lst_quants = [0.1,0.5,0.9]):
    #ds_mrms_biascorrected_filled, ds_mrms, ds_ref_proceeding, ds_correction_to_mrms, ds_ref_where_mrms_is_0_and_ref_is_not,lst_quants
    # quantiles of correction factor
    lst_new_data_arrays = ["mean_daily_correction_factor", "max_daily_correction_factor",
                            "mrms_nonbiascorrected_daily_totals_mm", "mrms_biascorrected_daily_totals_mm",
                                "total_ref_fillvalues_mm", "frac_of_tot_biascrctd_rain_from_ref_fill",
                                "hours_of_ref_fillvalues", "ref_daily_totals_mm", "mrms_biascorrected_minus_ref_mm",
                                    "mrms_nonbiascorrected_minus_ref_mm"]
    # where the non-bias corrected mrms dataset is zero, assign np.nan to the bias correction factor
    total_mb_mrms, dic_chunks_mrms = estimate_chunk_memory(ds_mrms)
    ds_correction_to_mrms = xr.where((ds_mrms == 0), x = np.nan, y = ds_correction_to_mrms).chunk(dic_chunks_mrms)
    # compute quantiles of bias correction factor
    da_quant = ds_correction_to_mrms.chunk(dict(time=-1, latitude = "auto", longitude = "auto")).rainrate.quantile(q=lst_quants, dim = "time", skipna = True)
    ds_mrms_biascorrected_filled["correction_factor_quantile"] = da_quant
    #
    ## write temporary file
    # bm_time = time.time()
    # tmp_zarr = f"{fldr_scratch_zarr}{in_date}_ds_mrms_biascorrected_filled2.zarr"
    # lst_tmp_files_to_delete.append(tmp_zarr)
    # # gc.collect()
    # ds_mrms_biascorrected_filled.chunk(dict(time = "auto", latitude = "auto", longitude = "auto")).to_zarr(tmp_zarr, mode = "w",
    #                                                                                                                     encoding = define_zarr_compression(ds_mrms_biascorrected_filled))
    # ds_mrms_biascorrected_filled = xr.open_zarr(store=tmp_zarr).chunk(dict(time = "auto", latitude = "auto", longitude = "auto"))
    # # gc.collect()
    # print("exported _ds_mrms_biascorrected_filled2 to zarr")
    # print(f"Time to export (min): {((time.time() - bm_time)/60):.2f} | total script runtime (min): {((time.time() - start_time)/60):.2f}")
    #
    #
    # mean correction factor
    ds_correction_daily_mean = ds_correction_to_mrms.mean("time", skipna = True)
    ds_mrms_biascorrected_filled["mean_daily_correction_factor"] = ds_correction_daily_mean.rainrate
    # min correction factor
    ds_correction_daily_min = ds_correction_to_mrms.min("time", skipna = True)
    ds_mrms_biascorrected_filled["min_daily_correction_factor"] = ds_correction_daily_min.rainrate
    # max correction factor
    ds_correction_daily_max = ds_correction_to_mrms.max("time", skipna = True)
    ds_mrms_biascorrected_filled["max_daily_correction_factor"] = ds_correction_daily_max.rainrate
    #
    # # gc.collect()
    # tmp_bias_crctd_fld = f"{fldr_scratch_zarr}{in_date}_bias_crctd_fld3.zarr"
    # lst_tmp_files_to_delete.append(tmp_bias_crctd_fld)
    # ds_mrms_biascorrected_filled.chunk(dict(time = "auto", latitude = "auto", longitude = "auto")).to_zarr(tmp_bias_crctd_fld, mode = "w")
    # ds_mrms_biascorrected_filled = xr.open_zarr(store=tmp_bias_crctd_fld).chunk(dict(time = "auto", latitude = "auto", longitude = "auto"))
    # print("exported scratch zarr with suffix _bias_crctd_fld3.zarr")
    #
    # computing daily total uncorrected mrms
    ds_mrms_biascorrected_filled["mrms_nonbiascorrected_daily_totals_mm"] = ds_mrms.rainrate.mean("time")*24
    # computing daily total corrected mrms
    ds_mrms_biascorrected_filled["mrms_biascorrected_daily_totals_mm"] = ds_mrms_biascorrected_filled.rainrate.mean("time")*24
    # total ref fill values for the day
    ds_ref_fillvals_daily_tot = ds_ref_where_mrms_is_0_and_ref_is_not.rainrate.mean("time")*24 # mean mm/hr per day times 24 hours in day
    #
    # bm_time = time.time()
    # # gc.collect()
    # tmp_bias_crctd_fld = f"{fldr_scratch_zarr}{in_date}_bias_crctd_fld4.zarr"
    # lst_tmp_files_to_delete.append(tmp_bias_crctd_fld)
    # ds_mrms_biascorrected_filled.chunk(dict(time = "auto", latitude = "auto", longitude = "auto")).to_zarr(tmp_bias_crctd_fld, mode = "w",
    #                                                                                                        encoding = define_zarr_compression(ds_mrms_biascorrected_filled))
    # ds_mrms_biascorrected_filled = xr.open_zarr(store=tmp_bias_crctd_fld).chunk(dict(time = "auto", latitude = "auto", longitude = "auto"))
    # # gc.collect()
    # print("exported scratch zarr with suffix _bias_crctd_fld4.zarr")
    # print(f"Time to export (min): {((time.time() - bm_time)/60):.2f} | total script runtime (min): {((time.time() - start_time)/60):.2f}")
    #
    ds_mrms_biascorrected_filled["total_ref_fillvalues_mm"] = ds_ref_fillvals_daily_tot
    #
    # # gc.collect()
    # tmp_bias_crctd_fld = f"{fldr_scratch_zarr}{in_date}_bias_crctd_fld5.zarr"
    # lst_tmp_files_to_delete.append(tmp_bias_crctd_fld)
    # ds_mrms_biascorrected_filled.chunk(dict(time = "auto", latitude = "auto", longitude = "auto")).to_zarr(tmp_bias_crctd_fld, mode = "w", encoding = define_zarr_compression(ds_mrms_biascorrected_filled))
    # ds_mrms_biascorrected_filled = xr.open_zarr(store=tmp_bias_crctd_fld).chunk(dict(time = "auto", latitude = "auto", longitude = "auto"))
    # print("exported scratch zarr with suffix _bias_crctd_fld5.zarr")
    #
    # fraction of daily total rain from ref fillvals
    ds_frac_tot_rain_from_fillvals = ds_mrms_biascorrected_filled["total_ref_fillvalues_mm"]/ds_mrms_biascorrected_filled["mrms_biascorrected_daily_totals_mm"]
    ds_mrms_biascorrected_filled["frac_of_tot_biascrctd_rain_from_ref_fill"] = ds_frac_tot_rain_from_fillvals
    #
    # bm_time = time.time()
    # # gc.collect()
    # tmp_bias_crctd_fld = f"{fldr_scratch_zarr}{in_date}_bias_crctd_fld6.zarr"
    # lst_tmp_files_to_delete.append(tmp_bias_crctd_fld)
    # ds_mrms_biascorrected_filled.chunk(dict(time = "auto", latitude = "auto", longitude = "auto")).to_zarr(tmp_bias_crctd_fld, mode = "w",
    #                                                                                                        encoding = define_zarr_compression(ds_mrms_biascorrected_filled))
    # ds_mrms_biascorrected_filled = xr.open_zarr(store=tmp_bias_crctd_fld).chunk(dict(time = "auto", latitude = "auto", longitude = "auto"))
    # print("exported scratch zarr with suffix _bias_crctd_fld6.zarr")
    # print(f"Time to export (min): {((time.time() - bm_time)/60):.2f} | total script runtime (min): {((time.time() - start_time)/60):.2f}")
    # gc.collect()
    #
    # total time ref fill values were used
    n_tsteps_of_ref_fill = xr.where(ds_ref_where_mrms_is_0_and_ref_is_not.rainrate>0, 1, 0).sum(dim='time')
    s_times = pd.Series(ds_ref_where_mrms_is_0_and_ref_is_not.time.values)
    tstep_hours = s_times.diff().mode()[0] / np.timedelta64(1, 'h')
    hrs_of_ref_fill = n_tsteps_of_ref_fill*tstep_hours
    hrs_of_ref_fill = xr.where(hrs_of_ref_fill>0, hrs_of_ref_fill, 0)
    ds_mrms_biascorrected_filled["hours_of_ref_fillvalues"] = hrs_of_ref_fill
    # computing daily total ref in mrms coordinates
    ## spatially resample stage iv to mrms resolution
    ds_ref_proceeding_to_mrms= spatial_resampling(ds_ref_proceeding, ds_mrms, "latitude", "longitude")
    ds_mrms_biascorrected_filled["ref_daily_totals_mm"] = ds_ref_proceeding_to_mrms.rainrate.mean("time")*24
    # computing differences in daily totals between ref and mrms (bias and nonbias corrected)
    ds_dif_crctd = ds_mrms_biascorrected_filled["mrms_biascorrected_daily_totals_mm"] - ds_mrms_biascorrected_filled["ref_daily_totals_mm"]
    ds_dif_uncrctd = ds_mrms_biascorrected_filled["mrms_nonbiascorrected_daily_totals_mm"] - ds_mrms_biascorrected_filled["ref_daily_totals_mm"]
    ds_mrms_biascorrected_filled["mrms_biascorrected_minus_ref_mm"] = ds_dif_crctd
    ds_mrms_biascorrected_filled["mrms_nonbiascorrected_minus_ref_mm"] = ds_dif_uncrctd
    return ds_mrms_biascorrected_filled, lst_new_data_arrays, lst_tmp_files_to_delete
# xds_mrms_biascorrected_filled= bias_correct_and_fill_mrms(ds_mrms, ds_ref)


#%%
tmp_raw_mrms_zarr = fldr_scratch_zarr + fl_in_zarr.split("/")[-1].split(".zarr")[0] + "_raw.zarr"

dic_auto_chunk = {'time':'auto', 'latitude': "auto", 'longitude': "auto"}
dic_mrms_chunks = {'time':1, 'latitude': 3500, 'longitude': 3500}

lst_tmp_files_to_delete = []
# try:
ds_mrms = xr.open_dataset(fl_in_zarr, chunks = dic_mrms_chunks, engine = "zarr")




total_mb_mrms, dic_chunks_mrms = estimate_chunk_memory(ds_mrms, dic_mrms_chunks)
# print(ds_mrms)
print("MRMS chunk memory (mb) and chunks: {:.2f}, {}".format(total_mb_mrms, dic_chunks_mrms))

performance["filepath_mrms"] = fl_in_zarr
# create a single row dataset with netcdf attributes
columns = []
values = []
columns.append('filepath')
values.append(ds_mrms.encoding['source'])
for key in ds_mrms.attrs:
    if isinstance(ds_mrms.attrs[key], dict):
        for key2 in ds_mrms.attrs[key]:
            columns.append(key2)
            values.append(str(ds_mrms.attrs[key][key2]))
    if key == 'warnings':
        input_string = ds_mrms.attrs[key]
        result_dict = ast.literal_eval(input_string)
        for key2 in result_dict:
            columns.append(key2)
            values.append(str(result_dict[key2]))
    else:
        columns.append(key)
        values.append(str(ds_mrms.attrs[key]))
# select subset based on the extents of the transposition domain
# the +360 is to convert from degrees west to degrees east; the + or - 0.05 is to buffer the selction by 5 gridcells assuming 0.01 degree grid
if gdf_transdomain is not None:
    ds_mrms = clip_ds_to_transposition_domain(ds_mrms, gdf_transdomain)
# replace na and negative values with 0
ds_mrms = ds_mrms.fillna(0)
ds_mrms = ds_mrms.where(ds_mrms>=0, 0, drop=False) # if negative values are present, replace them with 0

# update mrms data to have centered gridcells if aorc is used for correction
if bias_correction_reference == "aorc":
    # reindex mrms gridcells to represent centers
    grid_spacing = ds_mrms.attrs["grid_spacing"]
    if ds_mrms.attrs["gridcell_feature_represented_by_coordinate"] == 'upper_left':
        ds_mrms["latitude"] = ds_mrms.latitude - grid_spacing/2 # shift DOWN to middel
        ds_mrms["longitude"] = ds_mrms.longitude + grid_spacing/2 # shift RIGHT to middle
        ds_mrms.attrs["gridcell_feature_represented_by_coordinate"] = "center"


bm_time = time.time()
lst_tmp_files_to_delete.append(tmp_raw_mrms_zarr)
ds_mrms.chunk(dict(dic_mrms_chunks)).to_zarr(tmp_raw_mrms_zarr, mode = "w", encoding = define_zarr_compression(ds_mrms))
ds_mrms = xr.open_dataset(tmp_raw_mrms_zarr, chunks = dic_mrms_chunks, engine = "zarr")
print(f"Time to export rechunked mrms dataset (min): {((time.time() - bm_time)/60):.2f} | total script runtime (min): {((time.time() - start_time)/60):.2f}")

# print("Loaded MRMS data and filled missing and negative values with 0")
if performance["data_available_for_bias_correction"]:
    if bias_correction_reference == "stageiv":
        performance["filepath_ref_data"] = f_nc_stageiv
        ds_ref = xr.open_dataset(f_nc_stageiv, chunks = dic_auto_chunk)
        ds_ref = process_dans_stageiv(ds_ref, ds_mrms)
        ds_ref['time'] = ds_ref.time - pd.Timedelta(1, "hours")

    if bias_correction_reference == "aorc":
        performance["filepath_ref_data"] = f_aorc_yr
        da_aorc_rainfall.name = "rainrate" # rename like mrms
        if da_aorc_rainfall.attrs["units"] != "kg/m^2": # verify units (these are equivalent to mm/hr)
            sys.exit(f'WARNING: AORC RAINFALL UNITS NOT RECOGNIZED da_aorc_rainfall.attrs["units"]={da_aorc_rainfall.attrs["units"]}')
        ds_ref = da_aorc_rainfall.to_dataset()
        if ds_ref.longitude.values.min() < 0:
            ds_ref["longitude"] = ds_ref.longitude + 360
        
    if gdf_transdomain is not None:
        ds_ref = clip_ds_to_transposition_domain(ds_ref, gdf_transdomain)
    # replace na and negative values with 0 (there shouldn't be any so this is just to make sure)
    ds_ref = ds_ref.fillna(0) 
    ds_ref = ds_ref.where(ds_ref>=0, 0, drop=False) # if negative values are present, replace them with 0
    #
    # write to zarr and re-load dataset
    bm_time = time.time()
    tmp_raw_ds_ref_zarr = fldr_scratch_zarr + in_date + "_ref_data_raw.zarr"
    lst_tmp_files_to_delete.append(tmp_raw_ds_ref_zarr)
    ds_ref.chunk(dic_mrms_chunks).to_zarr(tmp_raw_ds_ref_zarr, mode = "w", encoding = define_zarr_compression(ds_ref))
    ds_ref = xr.open_dataset(tmp_raw_ds_ref_zarr, chunks = dic_mrms_chunks, engine = "zarr")
    # print("Loaded Stage IV data and filled missing and negative values with 0")
    print(f"Wrote raw {bias_correction_reference} data to zarr: {((time.time() - bm_time)/60):.2f} | total script runtime (min): {((time.time() - start_time)/60):.2f}")

    ds_mrms_biascorrected_filled,ds_mrms_hourly_to_ref,ds_ref_proceeding,\
            ds_correction_to_mrms, ds_ref_where_mrms_is_0_and_ref_is_not,\
                    lst_tmp_files_to_delete = bias_correct_and_fill_mrms(ds_mrms, ds_ref, lst_tmp_files_to_delete)
    
    # ds_mrms_biascorrected_filled,ds_mrms_hourly_to_ref,ds_ref_proceeding, ds_correction_to_mrms, ds_ref_where_mrms_is_0_and_ref_is_not, lst_tmp_files_to_delete = xds_mrms_biascorrected_filled, xds_mrms_hourly_to_ref, ds_ref, xds_correction_to_mrms, xds_ref_where_mrms_is_0_and_ref_is_not, lst_tmp_files_to_delete
    # print("ran function bias_correct_and_fill_mrms")
    # ds_mrms_biascorrected_filled,ds_mrms_hourly_to_ref,ds_ref_proceeding,\
    #         ds_correction_to_mrms, ds_ref_where_mrms_is_0_and_ref_is_not = xds_mrms_biascorrected_filled, xds_mrms_hourly_to_ref, ds_ref, xds_correction_to_mrms, xds_ref_where_mrms_is_0_and_ref_is_not
    # # gc.collect()
    # tmp_bias_crctd_fld = f"{fldr_scratch_zarr}{in_date}_bias_crctd_fld2.zarr"
    # lst_tmp_files_to_delete.append(tmp_bias_crctd_fld)
    # ds_mrms_biascorrected_filled.chunk(dict(time = "auto", latitude = "auto", longitude = "auto")).to_zarr(tmp_bias_crctd_fld, mode = "w")
    # ds_mrms_biascorrected_filled = xr.open_zarr(store=tmp_bias_crctd_fld).chunk(dict(time = "auto", latitude = "auto", longitude = "auto"))
    # print("exported bias corrected and filled mrms dataset to zarr (again)")
    # tmp0_ds_biascorrected_filled_zarr = fldr_scratch_zarr + fl_in_nc.split("/")[-1].split(".nc")[0] + "_processed0.zarr"
    # lst_tmp_files_to_delete.append(tmp0_ds_biascorrected_filled_zarr)
    # print("exporting intermediate output........")
    # time_size = ds_mrms_biascorrected_filled.sizes['time']
    # ds_mrms_biascorrected_filled.chunk(dict(time = "auto", latitude = "auto", longitude = "auto")).to_zarr(tmp0_ds_biascorrected_filled_zarr, mode = "w")
    # print("exported temporary bias corrected dataset to zarr (first intermediate output)")
    # ds_mrms_biascorrected_filled = xr.open_zarr(store=tmp0_ds_biascorrected_filled_zarr).chunk(dict(time = "auto", latitude = "auto", longitude = "auto"))
    # print("loaded temporary bias corrected dataset from zarr to consolidate to targeted timestep (first intermediate output)")
    # print("ran function bias_correct_and_fill_mrms")
    ds_mrms_biascorrected_filled,lst_new_data_arrays, lst_tmp_files_to_delete = process_bias_corrected_dataset(ds_mrms_biascorrected_filled, ds_mrms, ds_ref_proceeding,
                                    ds_correction_to_mrms, ds_ref_where_mrms_is_0_and_ref_is_not, lst_tmp_files_to_delete,
                                    lst_quants)

    # print("ran function process_bias_corrected_dataset")
    performance["domainwide_totals_CORRECTED_mrms_mm"] = ds_mrms_biascorrected_filled.attrs["domainwide_totals_CORRECTED_mrms_mm"]
    performance["domainwide_totals_uncorrected_mrms_mm"] = ds_mrms_biascorrected_filled.attrs["domainwide_totals_uncorrected_mrms_mm"]
    performance["domainwide_totals_ref_mm"] = ds_mrms_biascorrected_filled.attrs["domainwide_totals_ref_mm"]
    # performance["domainwide_totals_CORRECTED_mrms_over_ref"] = ds_mrms_biascorrected_filled.attrs["domainwide_totals_CORRECTED_mrms_over_ref"]
    performance["domainwide_totals_uncorrected_mrms_over_ref"] = ds_mrms_biascorrected_filled.attrs["domainwide_totals_uncorrected_mrms_over_ref"]
    performance["correction_factor_upperbound"] = ds_mrms_biascorrected_filled.attrs["correction_factor_upperbound"]
    performance["correction_factor_lowerbound"] = ds_mrms_biascorrected_filled.attrs["correction_factor_lowerbound"]
    print("finished bias correcting and filling mrms dataset")
else:
    performance["loading_raw_data_errors"] = False
df_input_dataset_attributes = pd.DataFrame([values], columns=columns)


# tstep = ds.attrs["time_step"]
if performance["data_available_for_bias_correction"]:
    ds_to_export = ds_mrms_biascorrected_filled
    # bm_time = time.time()
    # tmp_ds_biascorrected_filled_zarr = fldr_scratch_zarr + fl_in_zarr.split("/")[-1].split(".zarr")[0] + "_processed.zarr"
    # lst_tmp_files_to_delete.append(tmp_ds_biascorrected_filled_zarr)
    # # gc.collect()
    # ds_to_export.chunk(dict(time = "auto", latitude = "auto", longitude = "auto")).to_zarr(tmp_ds_biascorrected_filled_zarr, mode = "w",
    #                                                                                         encoding = define_zarr_compression(ds_to_export))
    # print("exported temporary bias corrected dataset to zarr")
    # ds_to_export = xr.open_zarr(store=tmp_ds_biascorrected_filled_zarr).chunk(dict(time = "auto", latitude = "auto", longitude = "auto"))
    # gc.collect()
    # print("loaded temporary bias corrected dataset from zarr to consolidate to targeted timestep")
    # print(f"Time to export (min): {((time.time() - bm_time)/60):.2f} | total script runtime (min): {((time.time() - start_time)/60):.2f}")
    # add dimension so data can be indexed based on whether it was bias corrected
    ds_to_export = ds_to_export.assign_coords(bias_corrected=True)
    ds_to_export = ds_to_export.expand_dims("bias_corrected")
else:
    ds_to_export = ds_mrms
    ds_to_export = ds_to_export.assign_coords(bias_corrected=False)
    ds_to_export = ds_to_export.expand_dims("bias_corrected")
# verify the full day has coverage
tstep_min = pd.to_timedelta(ds_mrms.attrs["time_step_min"]).total_seconds() / 60
num_tsteps = ds_to_export.coords["time"].shape[0]
duration_h = num_tsteps * tstep_min / 60
performance["duration_h"] = duration_h
performance["problem_with_duration"] = False
if duration_h != 24:
    performance["problem_with_duration"] = True
performance["current_tstep_different_than_target"] = False
if tstep_min != target_tstep_min: # consolidate to target timestep
    performance["current_tstep_different_than_target"] = True
    # resampling
    # performance["problems_resampling"] = True
    t_idx_1min = pd.date_range(ds_to_export.time.values[0], periods = 24*60, freq='1min')
    da_1min = ds_to_export.rainrate.reindex(dict(time = t_idx_1min)).ffill(dim="time")
    da_target = da_1min.resample(time = "{}Min".format(target_tstep_min)).mean()
    del ds_to_export['rainrate']
    ds_to_export['time'] = da_target.time
    ds_to_export['rainrate'] = da_target
    # performance["problems_resampling"] = False
try:
    bm_time = time.time()
    tmp_ds_penultimate_zarr = fldr_scratch_zarr + fl_in_zarr.split("/")[-1].split(".zarr")[0] + "_penultimate.zarr"
    lst_tmp_files_to_delete.append(tmp_ds_penultimate_zarr)
    ds_to_export.chunk(dict(time = "auto", latitude = "auto", longitude = "auto")).to_zarr(tmp_ds_penultimate_zarr, mode = "w",
                                                                                            encoding = define_zarr_compression(ds_to_export) , consolidated=True)
    ds_to_export = xr.open_zarr(store=tmp_ds_penultimate_zarr).chunk(dict(time = "auto", latitude = "auto", longitude = "auto"))
    print(f"Time to export penultimate zarr (min): {((time.time() - bm_time)/60):.2f} | total script runtime (min): {((time.time() - start_time)/60):.2f}")
    bm_time = time.time()
    # gc.collect()
    # chunk based on rainrate
    size, chnk_dic = estimate_chunk_memory(ds_to_export["rainrate"], input_chunk_sizes=final_chunking_dict)
    print(f"exporting to zarr with chunking {chnk_dic}")
    ds_to_export.chunk(chnk_dic).to_zarr(fl_out_zarr, mode="w", encoding=define_zarr_compression(ds_to_export), consolidated=True)
    # gc.collect()
    print(f"time to export zarr (min): {((time.time() - bm_time)/60):.2f} | total script runtime (min): {((time.time() - start_time)/60):.2f}")
    performance["to_zarr_errors"] = "None"
    performance["problem_exporting_zarr"] = False
    if final_output_type == "zarr":
        print(f"Wrote {fl_out_zarr}")
except Exception as e:
    performance["to_zarr_errors"]  = e
    performance["problem_exporting_zarr"] = True

if (final_output_type == "nc") and (performance["problem_exporting_zarr"] == False):
    try:
        # print("exporting to netcdf....")
        bm_time = time.time()
        fl_out_nc = fl_out_zarr.replace("zarr", "nc")
        Path(fl_out_nc).parent.mkdir(parents=True, exist_ok=True)
        ds_to_export = xr.open_zarr(store=fl_out_zarr).chunk(dic_mrms_chunks)
        encoding={var: {"zstd": True, "complevel": 5} for var in ds_to_export.data_vars}
        ds_to_export.to_netcdf(fl_out_nc, encoding = encoding, engine = "h5netcdf")
        print(f"time to export netcdf (min): {((time.time() - bm_time)/60):.2f} | total script runtime (min): {((time.time() - start_time)/60):.2f}")
        # remove the zarr file
        shutil.rmtree(fl_out_zarr)
        performance["to_nc_errors"] = "None"
        performance["problem_exporting_nc"] = False
        print(f"Wrote {fl_out_nc}")
    except Exception as e:
        performance["to_nc_errors"] = "None"
        performance["problem_exporting_nc"] = False

# delete intermediate inputs:
for f in lst_tmp_files_to_delete:
    shutil.rmtree(f)

# export performance dictionary to a csv
time_elapsed_min = round((time.time() - start_time) / 60, 2)
performance["time_elapsed_min"] = time_elapsed_min
df = pd.DataFrame(performance, index = [1])
df.to_csv(fl_out_csv)
print("wrote file: {}".format(fl_out_csv))
df_input_dataset_attributes.to_csv(fl_out_csv_qaqc)
print("wrote file: {}".format(fl_out_csv_qaqc))
print(f"script finished. Elapsed time (min): {performance['time_elapsed_min']}")
# %%
