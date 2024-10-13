#%% Import libraries
import time
import shutil
import xarray as xr
import geopandas as gp
import pandas as pd
import zarr
import sys
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



target_tstep = 2

# chnk_sz = "100MB"

performance = {}
#%% work

# in_date = "20210720" # "20210719" corresponds to slurm task array 200
# fldr_zarr_fullres_daily = "/project/quinnlab/dcl3nd/norfolk/highres-radar-rainfall-processing/data/mrms_zarr_preciprate_fullres_dailyfiles/"
# fldr_zarr_fullres_daily_constant_tstep = "/project/quinnlab/dcl3nd/norfolk/highres-radar-rainfall-processing/data/mrms_zarr_preciprate_fullres_dailyfiles_constant_tstep/"
# fldr_scratch_zarr = "/project/quinnlab/dcl3nd/norfolk/highres-radar-rainfall-processing/data/_scratch/zarrs/"
# fldr_scratch_csv = "/project/quinnlab/dcl3nd/norfolk/highres-radar-rainfall-processing/data/_scratch/csv/"
# fldr_nc_stageiv = "/project/quinnlab/dcl3nd/norfolk/stormy/data/climate/StageIV_rainfall/"
# # f_shp_sst_transom = "/project/quinnlab/dcl3nd/stormy/stochastic_storm_transposition/norfolk/transposition_domain/norfolk_trans_dom_4326.shp"
# f_shp_sst_transom = None
#%% end work



# folders (with proceeding fwd slash)
in_date = str(sys.argv[1]) # YYYYMMDD
fldr_zarr_fullres_daily = str(sys.argv[2]) # ${assar_dirs[out_fullres_dailyfiles]} # "/project/quinnlab/dcl3nd/norfolk/highres-radar-rainfall-processing/data/mrms_zarr_preciprate_fullres_dailyfiles/"
fldr_zarr_fullres_daily_constant_tstep = str(sys.argv[3]) # ${assar_dirs[out_fullres_dailyfiles_consolidated]} # "/project/quinnlab/dcl3nd/highres-radar-rainfall-processing/mrms_zarr_preciprate_fullres_dailyfiles_constant_tstep/"
fldr_scratch_zarr = str(sys.argv[4]) # ${assar_dirs[scratch_zarrs]} # "/project/quinnlab/dcl3nd/norfolk/highres-radar-rainfall-processing/data/_scratch/zarrs/"
fldr_scratch_csv = str(sys.argv[5]) # ${assar_dirs[scratch_zarrs]} # "/project/quinnlab/dcl3nd/norfolk/highres-radar-rainfall-processing/data/_scratch/csv/"
fldr_nc_stageiv = str(sys.argv[6])
try:
    f_shp_sst_transom = str(sys.argv[7]) # ${assar_dirs[shp_transposition_domain]} # "/project/quinnlab/dcl3nd/norfolk/stormy/stochastic_storm_transposition/norfolk/transposition_domain/norfolk_trans_dom_4326.shp"
except:
    f_shp_sst_transom = None
    print("Bias correcting entire MRMS dataset since no transposition domain shapefile has been supplied........")
performance["date"] = in_date

# f_out_export_perf = fldr_scratch_zarr + "_export_stats_{}.csv".format(in_date)
#%%  
fl_in_zarr = fldr_zarr_fullres_daily +"{}.zarr".format(in_date)
fl_out_zarr = fldr_zarr_fullres_daily_constant_tstep +"{}.zarr".format(in_date)
# fl_out_zarr = fldr_scratch_zarr +"{}.zarr".format(in_date)
fl_out_csv = fldr_scratch_csv +"da2_resampling_{}.csv".format(in_date)
fl_out_csv_qaqc = fldr_scratch_csv +"qaqc_of_daily_fullres_data_{}.csv".format(in_date)
# find related stageIV rainfall file
f_nc_stageiv = glob(fldr_nc_stageiv + "{}/*{}*".format(in_date[0:4],in_date))
if len(f_nc_stageiv)>0:
    f_nc_stageiv = f_nc_stageiv[0]
    stageiv_data_available_for_bias_correction = True
else:
    stageiv_data_available_for_bias_correction = False
    print("No stage iv data for this date. No bias correction being performed....")

if f_shp_sst_transom is not None:
    gdf_transdomain = gp.read_file(f_shp_sst_transom)
else: 
    gdf_transdomain = None

def define_zarr_compression(ds, clevel=5):
    encoding = {}
    for da_name in ds.data_vars:
        encoding[da_name] = {"compressor": zarr.Blosc(cname="zlib", clevel=clevel, shuffle=zarr.Blosc.SHUFFLE)}
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


def bias_correct_and_fill_mrms(ds_mrms, ds_stageiv, lst_tmp_files_to_delete, 
                               crxn_upper_bound = crxn_upper_bound, crxn_lower_bound = crxn_lower_bound,
                                 verbose = False):
    # convert stage iv to proceeding time interval
    # ds_stageiv_og = xr.open_dataset(tmp_raw_stage_iv_zarr, chunks = dic_chunks, engine = "zarr")
    # ds_stageiv = ds_stageiv_og.copy()
    ds_stageiv['time'] = ds_stageiv.time - pd.Timedelta(1, "hours")
    # create hourly version of mrms data
    ds_mrms_hourly = ds_mrms.resample(time = "H").mean() # convert to hourly timestep
    # spatially resample MRMS data to match stage IV resolution
    xds_mrms_hourly_to_stageiv= spatial_resampling(ds_mrms_hourly, ds_stageiv, "latitude", "longitude")
    #
    ## write the bias correction dataset to a temporary file
    tmp_zarr = f"{fldr_scratch_zarr}{in_date}_xds_mrms_hourly_to_stageiv.zarr"
    lst_tmp_files_to_delete.append(tmp_zarr)
    gc.collect()
    bm_time = time.time()
    xds_mrms_hourly_to_stageiv.chunk(dict(time = "auto", latitude = "auto", longitude = "auto")).to_zarr(tmp_zarr, mode = "w", encoding = define_zarr_compression(xds_mrms_hourly_to_stageiv))
    xds_mrms_hourly_to_stageiv = xr.open_zarr(store=tmp_zarr).chunk(dict(time = "auto", latitude = "auto", longitude = "auto"))
    gc.collect()
    print("exported xds_mrms_hourly_to_stageiv to zarr")
    print(f"Time to export (min): {((time.time() - bm_time)/60):.2f} | total script runtime (min): {((time.time() - start_time)/60):.2f}")
    #
    #
    # compute correction factor
    xds_mrms_hourly_correction_factor_st4res = ds_stageiv/xds_mrms_hourly_to_stageiv
    # if mrms is 0, assign a 1 correction factor (the case where stage iv is non zero is taken care of below) (otherwise it will be infinity)
    xds_mrms_hourly_correction_factor_st4res = xr.where((xds_mrms_hourly_to_stageiv <= 0),
                                                        x = 1, y = xds_mrms_hourly_correction_factor_st4res)
    # to mitigate outliers, enforce correction factor bounds
    ## where stage iv is 0 or negative, assign a value of 1 (so no correction, assuming stage iv is missing in these locations)
    ### this should only affect zeros since I'm covnerting negative values to 0
    xds_mrms_hourly_correction_factor_st4res = xr.where((ds_stageiv <= 0),
                                                        x = 1, y = xds_mrms_hourly_correction_factor_st4res)
    ## where upper bound is exceeded, assign upper bound
    xds_mrms_hourly_correction_factor_st4res = xr.where((xds_mrms_hourly_correction_factor_st4res > crxn_upper_bound),
                                                        x = crxn_upper_bound, y = xds_mrms_hourly_correction_factor_st4res)
    ## where lower bound is exceeded, assign lower bound
    xds_mrms_hourly_correction_factor_st4res = xr.where((xds_mrms_hourly_correction_factor_st4res < crxn_lower_bound),
                                                        x = crxn_lower_bound, y = xds_mrms_hourly_correction_factor_st4res)
    # upsample correction factor in space to MRMS resolution
    xds_mrms_hourly_correction_factor_fulres= spatial_resampling(xds_mrms_hourly_correction_factor_st4res, ds_mrms_hourly, "latitude", "longitude")
    #
    ## write the bias correction dataset to a temporary file
    tmp_zarr = f"{fldr_scratch_zarr}{in_date}_xds_mrms_hourly_correction_factor_fulres.zarr"
    lst_tmp_files_to_delete.append(tmp_zarr)
    bm_time = time.time()
    gc.collect()
    xds_mrms_hourly_correction_factor_fulres.chunk(dict(time = "auto", latitude = "auto", longitude = "auto")).to_zarr(tmp_zarr, mode = "w",
                                                                                                                        encoding = define_zarr_compression(xds_mrms_hourly_correction_factor_fulres))
    xds_mrms_hourly_correction_factor_fulres = xr.open_zarr(store=tmp_zarr).chunk(dict(time = "auto", latitude = "auto", longitude = "auto"))
    gc.collect()
    print("exported xds_mrms_hourly_correction_factor_fulres to zarr")
    print(f"Time to export (min): {((time.time() - bm_time)/60):.2f} | total script runtime (min): {((time.time() - start_time)/60):.2f}")
    #
    #
    # if stageiv is non zero and mrms is zero, fill with stage iv precip intensities
    ### create stage iv data with same res as mrms
    xds_stageiv_to_mrms= spatial_resampling(ds_stageiv, ds_mrms_hourly, "latitude", "longitude")
    #
    ## write temporary file
    bm_time = time.time()
    tmp_zarr = f"{fldr_scratch_zarr}{in_date}_xds_stageiv_to_mrms.zarr"
    lst_tmp_files_to_delete.append(tmp_zarr)
    gc.collect()
    xds_stageiv_to_mrms.chunk(dict(time = "auto", latitude = "auto", longitude = "auto")).to_zarr(tmp_zarr, mode = "w",
                                                                                                                        encoding = define_zarr_compression(xds_stageiv_to_mrms))
    xds_stageiv_to_mrms = xr.open_zarr(store=tmp_zarr).chunk(dict(time = "auto", latitude = "auto", longitude = "auto"))
    gc.collect()
    print("exported xds_stageiv_to_mrms to zarr")
    print(f"Time to export (min): {((time.time() - bm_time)/60):.2f} | total script runtime (min): {((time.time() - start_time)/60):.2f}")
    #
    #
    ### create dataset with stageiv rainfall intensities at the indieces where condition is true and 0 everywhere else
    xds_stage_iv_where_mrms_is_0_and_stageiv_is_not_hourly = xr.where((xds_stageiv_to_mrms > 0) & (ds_mrms_hourly == 0),
                                                            x = xds_stageiv_to_mrms, y = 0)
    #### upsample to the full MRMS resolution
    ##### forward fill missing values since it is a proceeding dataset
    xds_stage_iv_where_mrms_is_0_and_stageiv_is_not = xds_stage_iv_where_mrms_is_0_and_stageiv_is_not_hourly.reindex(dict(time = ds_mrms.time)).ffill(dim="time").chunk(dict(latitude = "auto", longitude = "auto"))
    print("exporting xds_stage_iv_where_mrms_is_0_and_stageiv_is_not (this is around a common kill point)")
    bm_time = time.time()
    tmp_zarr= f"{fldr_scratch_zarr}{in_date}_xds_stage_iv_where_mrms_is_0_and_stageiv_is_not.zarr"
    lst_tmp_files_to_delete.append(tmp_zarr)
    gc.collect()
    xds_stage_iv_where_mrms_is_0_and_stageiv_is_not.chunk(dict(time = -1, latitude = "auto", longitude = "auto")).to_zarr(tmp_zarr, mode = "w", encoding = define_zarr_compression(xds_stage_iv_where_mrms_is_0_and_stageiv_is_not))
    xds_stage_iv_where_mrms_is_0_and_stageiv_is_not = xr.open_zarr(store=tmp_zarr).chunk(dict(time = -1, latitude = "auto", longitude = "auto"))
    gc.collect()
    print("exported xds_stage_iv_where_mrms_is_0_and_stageiv_is_not to zarr")
    print(f"Time to export (min): {((time.time() - bm_time)/60):.2f} | total script runtime (min): {((time.time() - start_time)/60):.2f}")
    
    ### upsample bias correction to full res data
    xds_correction_to_mrms = xds_mrms_hourly_correction_factor_fulres.chunk(dict(time = -1, latitude = "auto", longitude = "auto")).reindex(dict(time = ds_mrms.time)).ffill(dim="time")
    print("exporting xds_correction_to_mrms (this is around a common kill point)")
    bm_time = time.time()
    tmp_zarr= f"{fldr_scratch_zarr}{in_date}_xds_correction_to_mrms.zarr"
    lst_tmp_files_to_delete.append(tmp_zarr)
    gc.collect()
    xds_correction_to_mrms.chunk(dict(time = -1, latitude = "auto", longitude = "auto")).to_zarr(tmp_zarr, mode = "w", encoding = define_zarr_compression(xds_correction_to_mrms))
    xds_correction_to_mrms = xr.open_zarr(store=tmp_zarr).chunk(dict(time = -1, latitude = "auto", longitude = "auto"))
    gc.collect()
    print("exported xds_correction_to_mrms to zarr")
    print(f"Time to export (min): {((time.time() - bm_time)/60):.2f} | total script runtime (min): {((time.time() - start_time)/60):.2f}")
    #
    #
    # write the bias correction dataset to a temporary file
    tmp_bias_correction_factor = f"{fldr_scratch_zarr}{in_date}_bias_crxn_factor.zarr"
    lst_tmp_files_to_delete.append(tmp_bias_correction_factor)
    gc.collect()
    # time_before_export = pd.Series(xds_correction_to_mrms.time.values)
    print("attempting to export xds_correction_to_mrms (common script kill point)")
    encoding = define_zarr_compression(xds_correction_to_mrms)
    encoding['time'] = {k: ds_mrms.time.encoding[k] for k in ['units', 'calendar', 'dtype'] if k in ds_mrms.time.encoding}
    xds_correction_to_mrms.chunk(dict(time = -1, latitude = "10MB", longitude = -1)).to_zarr(tmp_bias_correction_factor, mode = "w", encoding = encoding)
    xds_correction_to_mrms = xr.open_zarr(store=tmp_bias_correction_factor).chunk(dict(time = -1, latitude = "auto", longitude = "auto"))
    # time_after_export = pd.Series(xds_correction_to_mrms.time.values)
    print("exported bias correction factor to zarr")
    print(f"Time to export (min): {((time.time() - bm_time)/60):.2f} | total script runtime (min): {((time.time() - start_time)/60):.2f}")
    gc.collect()
    #
    ### apply correction factor
    xds_mrms_biascorrected = (ds_mrms * xds_correction_to_mrms)
    #
    bm_time = time.time()
    tmp_bias_crctd = f"{fldr_scratch_zarr}{in_date}_bias_crctd.zarr"
    lst_tmp_files_to_delete.append(tmp_bias_crctd)
    gc.collect()
    xds_mrms_biascorrected.chunk(dict(time = -1, latitude = "auto", longitude = "auto")).to_zarr(tmp_bias_crctd, mode = "w", encoding = define_zarr_compression(xds_mrms_biascorrected))
    xds_mrms_biascorrected = xr.open_zarr(store=tmp_bias_crctd).chunk(dict(time = -1, latitude = "auto", longitude = "auto"))
    gc.collect()
    print("exported xds_mrms_biascorrected to zarr")
    print(f"Time to export (min): {((time.time() - bm_time)/60):.2f} | total script runtime (min): {((time.time() - start_time)/60):.2f}")
    #
    ### fill in with stageIV data where mrms data is missing
    xds_mrms_biascorrected_filled = xds_mrms_biascorrected + xds_stage_iv_where_mrms_is_0_and_stageiv_is_not
    #
    bm_time = time.time()
    tmp_bias_crctd_fld = f"{fldr_scratch_zarr}{in_date}_bias_crctd_fld.zarr"
    lst_tmp_files_to_delete.append(tmp_bias_crctd_fld)
    gc.collect()
    xds_mrms_biascorrected_filled.chunk(dict(time = "auto", latitude = "auto", longitude = "auto")).to_zarr(tmp_bias_crctd_fld, mode = "w", encoding = define_zarr_compression(xds_mrms_biascorrected_filled))
    xds_mrms_biascorrected_filled = xr.open_zarr(store=tmp_bias_crctd_fld).chunk(dict(time = "auto", latitude = "auto", longitude = "auto"))
    gc.collect()
    print("exported bias corrected and filled mrms dataset to zarr")
    print(f"Time to export (min): {((time.time() - bm_time)/60):.2f} | total script runtime (min): {((time.time() - start_time)/60):.2f}")
    #
    ### keep original mrms data
    # xds_mrms_biascorrected_filled = xds_mrms_biascorrected_filled.assign(rainrate_uncorrected = ds_mrms.rainrate)
    ### include bias correction ds
    # xds_mrms_biascorrected_filled = xds_mrms_biascorrected_filled.assign(mrms_bias_correction_multiplier = xds_correction_to_mrms.rainrate)
    ### include stageiv fill values
    # xds_mrms_biascorrected_filled = xds_mrms_biascorrected_filled.assign(stageiv_fillvals_where_mrms_is_0_and_stageiv_is_not = xds_stage_iv_where_mrms_is_0_and_stageiv_is_not.rainrate)
    # check
    ## compare stage iv with bias corrected total rainfall for the whole day 
    ### computing daily totals by finding the average daily intensity in mm per hour with the .mean("time") function and then multiplying by 24 hours
    #### comparing using mrms resolution
    tot_rain_mrms_corrected = compute_total_rainfall_over_domain(xds_mrms_biascorrected_filled.chunk(dict(time = "auto", latitude = "auto", longitude = "auto")))
    tot_rain_stageiv = compute_total_rainfall_over_domain(ds_stageiv)
    tot_rain_mrms_uncorrected = compute_total_rainfall_over_domain(ds_mrms_hourly)
    ## old method
    # tot_rain_mrms_corrected = (xds_mrms_biascorrected_filled.mean("time")*24).rainrate.sum().values
    # tot_rain_stageiv = (xds_stageiv_to_mrms.mean("time")*24).rainrate.sum().values
    # tot_rain_mrms_uncorrected = (ds_mrms.mean("time")*24).rainrate.sum().values
    # print("Fraction of domain-wide rainfall totals:")
    # print("Bias corrected MRMS data over stage iv data: {}".format(tot_rain_mrms_corrected/tot_rain_stageiv))
    # print("UNbias corrected MRMS data over stage iv data: {}".format(tot_rain_mrms_uncorrected/tot_rain_stageiv))
    #### comparing using stageiv resolution
    ### resample bias corrected mrms data to stage iv resolution to compare
    # xds_mrms_biascorrected_to_stageiv= spatial_resampling(xds_mrms_biascorrected_filled, ds_stageiv, "latitude", "longitude")
    # tot_rain_mrms_corrected = (xds_mrms_biascorrected_to_stageiv.mean("time")*24).rainrate.sum().values
    # tot_rain_mrms_uncorrected = (xds_mrms_biascorrected_to_stageiv.mean("time")*24).rainrate_uncorrected.sum().values
    # tot_rain_stageiv = (ds_stageiv.mean("time")*24).rainrate.sum().values
    # print("Fraction of domain-wide rainfall totals:")
    # print("Bias corrected MRMS data over stage iv data: {}".format(tot_rain_mrms_corrected/tot_rain_stageiv))
    # print("UNbias corrected MRMS data over stage iv data: {}".format(tot_rain_mrms_uncorrected/tot_rain_stageiv))
    if tot_rain_stageiv > 0:
        domainwide_totals_CORRECTED_mrms_over_stageiv = tot_rain_mrms_corrected/tot_rain_stageiv
        domainwide_totals_uncorrected_mrms_over_stageiv = tot_rain_mrms_uncorrected/tot_rain_stageiv
    else:
        domainwide_totals_CORRECTED_mrms_over_stageiv = domainwide_totals_uncorrected_mrms_over_stageiv = np.nan
    xds_mrms_biascorrected_filled.attrs["domainwide_totals_CORRECTED_mrms_mm"] = tot_rain_mrms_corrected
    xds_mrms_biascorrected_filled.attrs["domainwide_totals_stageiv_mm"] = tot_rain_stageiv
    xds_mrms_biascorrected_filled.attrs["domainwide_totals_uncorrected_mrms_mm"] = tot_rain_mrms_uncorrected
    xds_mrms_biascorrected_filled.attrs["domainwide_totals_CORRECTED_mrms_over_stageiv"] = domainwide_totals_CORRECTED_mrms_over_stageiv
    xds_mrms_biascorrected_filled.attrs["domainwide_totals_uncorrected_mrms_over_stageiv"] = domainwide_totals_uncorrected_mrms_over_stageiv
    xds_mrms_biascorrected_filled.attrs["correction_factor_upperbound"] = crxn_upper_bound
    xds_mrms_biascorrected_filled.attrs["correction_factor_lowerbound"] = crxn_lower_bound
    return xds_mrms_biascorrected_filled, xds_mrms_hourly_to_stageiv, ds_stageiv,\
            xds_correction_to_mrms, xds_stage_iv_where_mrms_is_0_and_stageiv_is_not, lst_tmp_files_to_delete

def process_bias_corrected_dataset(ds_mrms_biascorrected_filled, ds_mrms, ds_stageiv_proceeding, ds_correction_to_mrms,
                                    ds_stage_iv_where_mrms_is_0_and_stageiv_is_not, lst_tmp_files_to_delete,
                                      lst_quants = [0.1,0.5,0.9]):
    #ds_mrms_biascorrected_filled, ds_mrms, ds_stageiv_proceeding, ds_correction_to_mrms, ds_stage_iv_where_mrms_is_0_and_stageiv_is_not,lst_quants
    # quantiles of correction factor
    lst_new_data_arrays = ["mean_daily_correction_factor", "max_daily_correction_factor",
                            "mrms_nonbiascorrected_daily_totals_mm", "mrms_biascorrected_daily_totals_mm",
                              "total_stageiv_fillvalues_mm", "frac_of_tot_biascrctd_rain_from_stageiv_fill",
                                "hours_of_stageiv_fillvalues", "stageiv_daily_totals_mm", "mrms_biascorrected_minus_stageiv_mm",
                                  "mrms_nonbiascorrected_minus_stageiv_mm"]
    # where the non-bias corrected mrms dataset is zero, assign np.nan to the bias correction factor
    ds_correction_to_mrms = xr.where((ds_mrms == 0), x = np.nan, y = ds_correction_to_mrms).chunk(dict(time = -1, latitude = "auto", longitude = "auto"))
    # compute quantiles of bias correction factor
    da_quant = ds_correction_to_mrms.rainrate.quantile(q=lst_quants, dim = "time", skipna = True)
    ds_mrms_biascorrected_filled["correction_factor_quantile"] = da_quant
    #
    ## write temporary file
    bm_time = time.time()
    tmp_zarr = f"{fldr_scratch_zarr}{in_date}_ds_mrms_biascorrected_filled2.zarr"
    lst_tmp_files_to_delete.append(tmp_zarr)
    gc.collect()
    ds_mrms_biascorrected_filled.chunk(dict(time = "auto", latitude = "auto", longitude = "auto")).to_zarr(tmp_zarr, mode = "w",
                                                                                                                        encoding = define_zarr_compression(ds_mrms_biascorrected_filled))
    ds_mrms_biascorrected_filled = xr.open_zarr(store=tmp_zarr).chunk(dict(time = "auto", latitude = "auto", longitude = "auto"))
    gc.collect()
    print("exported _ds_mrms_biascorrected_filled2 to zarr")
    print(f"Time to export (min): {((time.time() - bm_time)/60):.2f} | total script runtime (min): {((time.time() - start_time)/60):.2f}")
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
    # gc.collect()
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
    # total stageiv fill values for the day
    ds_stageiv_fillvals_daily_tot = ds_stage_iv_where_mrms_is_0_and_stageiv_is_not.rainrate.mean("time")*24 # mean mm/hr per day times 24 hours in day
    #
    bm_time = time.time()
    gc.collect()
    tmp_bias_crctd_fld = f"{fldr_scratch_zarr}{in_date}_bias_crctd_fld4.zarr"
    lst_tmp_files_to_delete.append(tmp_bias_crctd_fld)
    ds_mrms_biascorrected_filled.chunk(dict(time = "auto", latitude = "auto", longitude = "auto")).to_zarr(tmp_bias_crctd_fld, mode = "w",
                                                                                                           encoding = define_zarr_compression(ds_mrms_biascorrected_filled))
    ds_mrms_biascorrected_filled = xr.open_zarr(store=tmp_bias_crctd_fld).chunk(dict(time = "auto", latitude = "auto", longitude = "auto"))
    gc.collect()
    print("exported scratch zarr with suffix _bias_crctd_fld4.zarr")
    print(f"Time to export (min): {((time.time() - bm_time)/60):.2f} | total script runtime (min): {((time.time() - start_time)/60):.2f}")
    #
    ds_mrms_biascorrected_filled["total_stageiv_fillvalues_mm"] = ds_stageiv_fillvals_daily_tot
    #
    # gc.collect()
    # tmp_bias_crctd_fld = f"{fldr_scratch_zarr}{in_date}_bias_crctd_fld5.zarr"
    # lst_tmp_files_to_delete.append(tmp_bias_crctd_fld)
    # ds_mrms_biascorrected_filled.chunk(dict(time = "auto", latitude = "auto", longitude = "auto")).to_zarr(tmp_bias_crctd_fld, mode = "w", encoding = define_zarr_compression(ds_mrms_biascorrected_filled))
    # ds_mrms_biascorrected_filled = xr.open_zarr(store=tmp_bias_crctd_fld).chunk(dict(time = "auto", latitude = "auto", longitude = "auto"))
    # print("exported scratch zarr with suffix _bias_crctd_fld5.zarr")
    #
    # fraction of daily total rain from stageiv fillvals
    ds_frac_tot_rain_from_fillvals = ds_mrms_biascorrected_filled["total_stageiv_fillvalues_mm"]/ds_mrms_biascorrected_filled["mrms_biascorrected_daily_totals_mm"]
    ds_mrms_biascorrected_filled["frac_of_tot_biascrctd_rain_from_stageiv_fill"] = ds_frac_tot_rain_from_fillvals
    #
    bm_time = time.time()
    gc.collect()
    tmp_bias_crctd_fld = f"{fldr_scratch_zarr}{in_date}_bias_crctd_fld6.zarr"
    lst_tmp_files_to_delete.append(tmp_bias_crctd_fld)
    ds_mrms_biascorrected_filled.chunk(dict(time = "auto", latitude = "auto", longitude = "auto")).to_zarr(tmp_bias_crctd_fld, mode = "w",
                                                                                                           encoding = define_zarr_compression(ds_mrms_biascorrected_filled))
    ds_mrms_biascorrected_filled = xr.open_zarr(store=tmp_bias_crctd_fld).chunk(dict(time = "auto", latitude = "auto", longitude = "auto"))
    print("exported scratch zarr with suffix _bias_crctd_fld6.zarr")
    print(f"Time to export (min): {((time.time() - bm_time)/60):.2f} | total script runtime (min): {((time.time() - start_time)/60):.2f}")
    gc.collect()
    #
    # total time stageiv fill values were used
    n_tsteps_of_stageiv_fill = xr.where(ds_stage_iv_where_mrms_is_0_and_stageiv_is_not.rainrate>0, 1, 0).sum(dim='time')
    s_times = pd.Series(ds_stage_iv_where_mrms_is_0_and_stageiv_is_not.time.values)
    tstep_hours = s_times.diff().mode()[0] / np.timedelta64(1, 'h')
    hrs_of_stageiv_fill = n_tsteps_of_stageiv_fill*tstep_hours
    hrs_of_stageiv_fill = xr.where(hrs_of_stageiv_fill>0, hrs_of_stageiv_fill, 0)
    ds_mrms_biascorrected_filled["hours_of_stageiv_fillvalues"] = hrs_of_stageiv_fill
    # computing daily total stageiv in mrms coordinates
    ## spatially resample stage iv to mrms resolution
    ds_stageiv_proceeding_to_mrms= spatial_resampling(ds_stageiv_proceeding, ds_mrms, "latitude", "longitude")
    ds_mrms_biascorrected_filled["stageiv_daily_totals_mm"] = ds_stageiv_proceeding_to_mrms.rainrate.mean("time")*24
    # computing differences in daily totals between stageiv and mrms (bias and nonbias corrected)
    ds_dif_crctd = ds_mrms_biascorrected_filled["mrms_biascorrected_daily_totals_mm"] - ds_mrms_biascorrected_filled["stageiv_daily_totals_mm"]
    ds_dif_uncrctd = ds_mrms_biascorrected_filled["mrms_nonbiascorrected_daily_totals_mm"] - ds_mrms_biascorrected_filled["stageiv_daily_totals_mm"]
    ds_mrms_biascorrected_filled["mrms_biascorrected_minus_stageiv_mm"] = ds_dif_crctd
    ds_mrms_biascorrected_filled["mrms_nonbiascorrected_minus_stageiv_mm"] = ds_dif_uncrctd
    return ds_mrms_biascorrected_filled, lst_new_data_arrays, lst_tmp_files_to_delete
# xds_mrms_biascorrected_filled= bias_correct_and_fill_mrms(ds_mrms, ds_stageiv)
#%%
tmp_raw_mrms_zarr = fldr_scratch_zarr + fl_in_zarr.split("/")[-1].split(".zarr")[0] + "_raw.zarr"
tmp_raw_stage_iv_zarr = fldr_scratch_zarr + f_nc_stageiv.split("/")[-1].split(".nc")[0] + "_raw.zarr"
# shutil.rmtree(tmp_raw_mrms_zarr)
# shutil.rmtree(tmp_raw_stage_iv_zarr)


dic_chunks = {'time':'auto', 'latitude': "auto", 'longitude': "auto"}
lst_tmp_files_to_delete = []
# try:
ds_mrms = xr.open_dataset(fl_in_zarr, chunks = dic_chunks, engine = "zarr")
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

bm_time = time.time()
lst_tmp_files_to_delete.append(tmp_raw_mrms_zarr)
gc.collect()
ds_mrms.chunk(dict(time = "auto", latitude = "auto", longitude = "auto")).to_zarr(tmp_raw_mrms_zarr, mode = "w", encoding = define_zarr_compression(ds_mrms))
ds_mrms = xr.open_dataset(tmp_raw_mrms_zarr, chunks = dic_chunks, engine = "zarr")
print("Exported mrms after filling missing with 0")
print(f"Time to export (min): {((time.time() - bm_time)/60):.2f} | total script runtime (min): {((time.time() - start_time)/60):.2f}")
gc.collect()

performance["stageiv_available_for_bias_correction"] = True

# print("Loaded MRMS data and filled missing and negative values with 0")
if stageiv_data_available_for_bias_correction:
    performance["filepath_stageiv"] = f_nc_stageiv
    ds_stageiv = xr.open_dataset(f_nc_stageiv, chunks = dic_chunks)
    ds_stageiv = process_dans_stageiv(ds_stageiv)
    if gdf_transdomain is not None:
        ds_stageiv = clip_ds_to_transposition_domain(ds_stageiv, gdf_transdomain)
    # replace na and negative values with 0 (there shouldn't be any so this is just to make sure)
    ds_stageiv = ds_stageiv.fillna(0) 
    ds_stageiv = ds_stageiv.where(ds_stageiv>=0, 0, drop=False) # if negative values are present, replace them with 0
    # write to zarr and re-load dataset
    bm_time = time.time()
    lst_tmp_files_to_delete.append(tmp_raw_stage_iv_zarr)
    gc.collect()
    ds_stageiv.chunk(dict(time = "auto", latitude = "auto", longitude = "auto")).to_zarr(tmp_raw_stage_iv_zarr, mode = "w", encoding = define_zarr_compression(ds_stageiv))
    ds_stageiv = xr.open_zarr(store=tmp_raw_stage_iv_zarr).chunk(dict(time = "auto", latitude = "auto", longitude = "auto"))
    ds_stageiv = xr.open_dataset(tmp_raw_stage_iv_zarr, chunks = dic_chunks, engine = "zarr")
    #
    print("Loaded Stage IV data and filled missing and negative values with 0")
    print(f"Time to export (min): {((time.time() - bm_time)/60):.2f} | total script runtime (min): {((time.time() - start_time)/60):.2f}")
    gc.collect()
    
    ds_mrms_biascorrected_filled,ds_mrms_hourly_to_stageiv,ds_stageiv_proceeding,\
            ds_correction_to_mrms, ds_stage_iv_where_mrms_is_0_and_stageiv_is_not,\
                    lst_tmp_files_to_delete = bias_correct_and_fill_mrms(ds_mrms, ds_stageiv, lst_tmp_files_to_delete)
    print("ran function bias_correct_and_fill_mrms")
    # ds_mrms_biascorrected_filled,ds_mrms_hourly_to_stageiv,ds_stageiv_proceeding,\
    #         ds_correction_to_mrms, ds_stage_iv_where_mrms_is_0_and_stageiv_is_not = xds_mrms_biascorrected_filled, xds_mrms_hourly_to_stageiv, ds_stageiv, xds_correction_to_mrms, xds_stage_iv_where_mrms_is_0_and_stageiv_is_not
    # gc.collect()
    # tmp_bias_crctd_fld = f"{fldr_scratch_zarr}{in_date}_bias_crctd_fld2.zarr"
    # lst_tmp_files_to_delete.append(tmp_bias_crctd_fld)
    # ds_mrms_biascorrected_filled.chunk(dict(time = "auto", latitude = "auto", longitude = "auto")).to_zarr(tmp_bias_crctd_fld, mode = "w")
    # ds_mrms_biascorrected_filled = xr.open_zarr(store=tmp_bias_crctd_fld).chunk(dict(time = "auto", latitude = "auto", longitude = "auto"))
    # print("exported bias corrected and filled mrms dataset to zarr (again)")
    # tmp0_ds_biascorrected_filled_zarr = fldr_scratch_zarr + fl_in_nc.split("/")[-1].split(".nc")[0] + "_processed0.zarr"
    # lst_tmp_files_to_delete.append(tmp0_ds_biascorrected_filled_zarr)
    # print("exporting intermediate output........")
    # time_size = ds_mrms_biascorrected_filled.sizes['time']
    # ds_mrms_biascorrected_filled.chunk(dict(time = time_size, latitude = -1, longitude = "auto")).to_zarr(tmp0_ds_biascorrected_filled_zarr, mode = "w")
    # print("exported temporary bias corrected dataset to zarr (first intermediate output)")
    # ds_mrms_biascorrected_filled = xr.open_zarr(store=tmp0_ds_biascorrected_filled_zarr).chunk(dict(time = "auto", latitude = "auto", longitude = "auto"))
    # print("loaded temporary bias corrected dataset from zarr to consolidate to targeted timestep (first intermediate output)")
    # print("ran function bias_correct_and_fill_mrms")
    ds_mrms_biascorrected_filled,lst_new_data_arrays, lst_tmp_files_to_delete = process_bias_corrected_dataset(ds_mrms_biascorrected_filled, ds_mrms, ds_stageiv_proceeding,
                                    ds_correction_to_mrms, ds_stage_iv_where_mrms_is_0_and_stageiv_is_not, lst_tmp_files_to_delete,
                                    lst_quants)
    print("ran function process_bias_corrected_dataset")
    performance["domainwide_totals_CORRECTED_mrms_mm"] = ds_mrms_biascorrected_filled.attrs["domainwide_totals_CORRECTED_mrms_mm"]
    performance["domainwide_totals_uncorrected_mrms_mm"] = ds_mrms_biascorrected_filled.attrs["domainwide_totals_uncorrected_mrms_mm"]
    performance["domainwide_totals_stageiv_mm"] = ds_mrms_biascorrected_filled.attrs["domainwide_totals_stageiv_mm"]
    # performance["domainwide_totals_CORRECTED_mrms_over_stageiv"] = ds_mrms_biascorrected_filled.attrs["domainwide_totals_CORRECTED_mrms_over_stageiv"]
    performance["domainwide_totals_uncorrected_mrms_over_stageiv"] = ds_mrms_biascorrected_filled.attrs["domainwide_totals_uncorrected_mrms_over_stageiv"]
    performance["correction_factor_upperbound"] = ds_mrms_biascorrected_filled.attrs["correction_factor_upperbound"]
    performance["correction_factor_lowerbound"] = ds_mrms_biascorrected_filled.attrs["correction_factor_lowerbound"]
    print("finished bias correcting and filling mrms dataset. Onto final processing steps...")
else:
    performance["loading_raw_data_errors"] = False
df_input_dataset_attributes = pd.DataFrame([values], columns=columns)


# tstep = ds.attrs["time_step"]
if stageiv_data_available_for_bias_correction:
    # ds_to_export = ds_mrms_biascorrected_filled
    bm_time = time.time()
    tmp_ds_biascorrected_filled_zarr = fldr_scratch_zarr + fl_in_zarr.split("/")[-1].split(".zarr")[0] + "_processed.zarr"
    lst_tmp_files_to_delete.append(tmp_ds_biascorrected_filled_zarr)
    gc.collect()
    ds_mrms_biascorrected_filled.chunk(dict(time = "auto", latitude = "auto", longitude = "auto")).to_zarr(tmp_ds_biascorrected_filled_zarr, mode = "w", encoding = define_zarr_compression(ds_mrms_biascorrected_filled))
    print("exported temporary bias corrected dataset to zarr")
    ds_to_export = xr.open_zarr(store=tmp_ds_biascorrected_filled_zarr).chunk(dict(time = "auto", latitude = "auto", longitude = "auto"))
    gc.collect()
    # add dimension so data can be indexed based on whether it was bias corrected
    ds_to_export = ds_to_export.assign_coords(bias_corrected=True)
    ds_to_export = ds_to_export.expand_dims("bias_corrected")
    print("loaded temporary bias corrected dataset from zarr to consolidate to targeted timestep")
    print(f"Time to export (min): {((time.time() - bm_time)/60):.2f} | total script runtime (min): {((time.time() - start_time)/60):.2f}")
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
if tstep_min != target_tstep: # consolidate to target timestep
    performance["current_tstep_different_than_target"] = True
    # resampling
    # performance["problems_resampling"] = True
    t_idx_1min = pd.date_range(ds_to_export.time.values[0], periods = 24*60, freq='1min')
    da_1min = ds_to_export.rainrate.reindex(dict(time = t_idx_1min)).ffill(dim="time")
    da_target = da_1min.resample(time = "{}Min".format(target_tstep)).mean()
    del ds_to_export['rainrate']
    ds_to_export['time'] = da_target.time
    ds_to_export['rainrate'] = da_target
    # performance["problems_resampling"] = False
performance["problem_exporting_zarr"] = False
performance["to_zarr_errors"] = "None"
try:
    print("exporting to zarr....")
    bm_time = time.time()
    gc.collect()
    ds_to_export.chunk(dict(time = "auto", latitude = "auto", longitude = "auto")).to_zarr(fl_out_zarr, mode="w", encoding=define_zarr_compression(ds_to_export))
    gc.collect()
    print(f"Time to export (min): {((time.time() - bm_time)/60):.2f} | total script runtime (min): {((time.time() - start_time)/60):.2f}")
    # ds_from_zarr = xr.open_zarr(store=fl_out_zarr).chunk(dict(time = "auto", latitude = "auto", longitude = "auto"))
except Exception as e:
    performance["to_zarr_errors"]  = e
    performance["problem_exporting_zarr"] = True

# delete intermediate inputs:
for f in lst_tmp_files_to_delete:
    shutil.rmtree(f)

# export performance dictionary to a csv
time_elapsed_min = round((time.time() - start_time) / 60, 2)
performance["time_elapsed_min"] = time_elapsed_min
df = pd.DataFrame(performance, index = [1])
df.to_csv(fl_out_csv)
print("wrote file: {}".format(fl_out_csv))
if performance["problem_loading_netcdf"] == False:
    df_input_dataset_attributes.to_csv(fl_out_csv_qaqc)
    print("wrote file: {}".format(fl_out_csv_qaqc))
print(f"script finished. Elapsed time (min): {performance['time_elapsed_min']}")
# %%
