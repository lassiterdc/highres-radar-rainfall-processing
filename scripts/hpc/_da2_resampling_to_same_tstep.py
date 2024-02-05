#%% Import libraries
import time
start_time = time.time()
import shutil
import xarray as xr
import geopandas as gp
import pandas as pd
import sys
import dask
dask.config.set(**{'array.slicing.split_large_chunks': False}) # to silence warnings of loading large slice into memory
dask.config.set(scheduler='synchronous') # this forces single threaded computations
from __utils import *
import ast
from glob import glob
from rasterio.enums import Resampling
import rioxarray

target_tstep = 5

chnk_sz = "1000MB"

performance = {}
#%% work

in_date = "20210719" # "20210719" corresponds to slurm task array 200
fldr_nc_fullres_daily = "/scratch/dcl3nd/highres-radar-rainfall-processing/data/mrms_nc_preciprate_fullres_dailyfiles/"
fldr_nc_fullres_daily_constant_tstep = "/scratch/dcl3nd/highres-radar-rainfall-processing/data/mrms_nc_preciprate_fullres_dailyfiles_constant_tstep/"
fldr_scratch_zarr = "/scratch/dcl3nd/highres-radar-rainfall-processing/data/_scratch/zarrs/"
fldr_scratch_csv = "/scratch/dcl3nd/highres-radar-rainfall-processing/data/_scratch/csv/"
f_shp_sst_transom = "/scratch/dcl3nd/stormy/stochastic_storm_transposition/norfolk/transposition_domain/norfolk_trans_dom_4326.shp"
fldr_nc_stageiv = "/scratch/dcl3nd/stormy/data/climate/StageIV_rainfall/"

#%% end work



# folders (with proceeding fwd slash)
in_date = str(sys.argv[1]) # YYYYMMDD
fldr_nc_fullres_daily = str(sys.argv[2]) # ${assar_dirs[out_fullres_dailyfiles]} # "/project/quinnlab/dcl3nd/norfolk/highres-radar-rainfall-processing/data/mrms_nc_preciprate_fullres_dailyfiles/"
fldr_nc_fullres_daily_constant_tstep = str(sys.argv[3]) # ${assar_dirs[out_fullres_dailyfiles_consolidated]} # "/scratch/dcl3nd/highres-radar-rainfall-processing/out_fullres_dailyfiles_consolidated/"
fldr_scratch_zarr = str(sys.argv[4]) # ${assar_dirs[scratch_zarrs]} # "/project/quinnlab/dcl3nd/norfolk/highres-radar-rainfall-processing/data/_scratch/zarrs/"
fldr_scratch_csv = str(sys.argv[5]) # ${assar_dirs[scratch_zarrs]} # "/project/quinnlab/dcl3nd/norfolk/highres-radar-rainfall-processing/data/_scratch/csv/"
f_shp_sst_transom = str(sys.argv[6]) # ${assar_dirs[shp_transposition_domain]} # "/project/quinnlab/dcl3nd/norfolk/stormy/stochastic_storm_transposition/norfolk/transposition_domain/norfolk_trans_dom_4326.shp"
fldr_nc_stageiv = str(sys.argv[7])

performance["date"] = in_date

# f_out_export_perf = fldr_scratch_zarr + "_export_stats_{}.csv".format(in_date)
#%% netcdf 
fl_in_nc = fldr_nc_fullres_daily +"{}.nc".format(in_date)
fl_out_nc = fldr_nc_fullres_daily_constant_tstep +"{}.nc".format(in_date)
fl_out_zarr = fldr_scratch_zarr +"{}.zarr".format(in_date)
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

performance["problem_loading_netcdf"] = False
performance["loading_netcdf_errors"]  = "None"

gdf_transdomain = gp.read_file(f_shp_sst_transom)

def clip_ds_to_transposition_domain(ds, gdf_transdomain):
    transdom_bounds = gdf_transdomain.bounds
    ds = ds.where((ds.latitude >= float(transdom_bounds.miny.iloc[0]-.05)) & (ds.latitude <= float(transdom_bounds.maxy.iloc[0]+.05)) & (ds.longitude >= float(transdom_bounds.minx.iloc[0]+360-.05)) & (ds.longitude <= float(transdom_bounds.maxx.iloc[0]+360+.05)), drop = True)
    return ds

def spatial_resampling(xds_to_resample, xds_target, lat_varname, lon_varname):
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
    return xds_to_resampled

def bias_correct_and_fill_mrms(ds_mrms, ds_stageiv_og, target_quant = 0.998):
    # convert stage iv to proceeding time interval
    ds_stageiv = ds_stageiv_og.copy()
    ds_stageiv['time'] = ds_stageiv.time - pd.Timedelta(1, "hours")
    # create hourly version of mrms data
    ds_mrms_hourly = ds_mrms.resample(time = "H").mean() # convert to hourly timestep
    xds_mrms_hourly_to_stageiv= spatial_resampling(ds_mrms_hourly, ds_stageiv, "latitude", "longitude")
    # compute correction factor
    xds_mrms_hourly_correction_factor_st4res = ds_stageiv/xds_mrms_hourly_to_stageiv
    # if mrms or stage iv are 0, assign a 0 correction factor (the case where stage iv is non zero is taken care of below)
    xds_mrms_hourly_correction_factor_st4res = xr.where((ds_stageiv == 0) | (xds_mrms_hourly_to_stageiv == 0),
                                                        x = 0, y = xds_mrms_hourly_correction_factor_st4res)
    # to mitigate crazy outliers:
    ## if the correction factor is greater than the 99th percentile correction factor, assign it the 99th percentile correction factor
    quant_correction = xds_mrms_hourly_correction_factor_st4res.to_dataframe().rainrate.quantile(target_quant)
    xds_mrms_hourly_correction_factor_st4res = xr.where((xds_mrms_hourly_correction_factor_st4res > quant_correction),
                                                        x = quant_correction, y = xds_mrms_hourly_correction_factor_st4res)
    xds_mrms_hourly_correction_factor_fulres= spatial_resampling(xds_mrms_hourly_correction_factor_st4res, ds_mrms_hourly, "latitude", "longitude")
    # replace the fillvalue with 0 (3.403e+38) - this corresponds to areas on the edges where MRMS falls outside the range of stage iv
    xds_mrms_hourly_correction_factor_fulres = xr.where(xds_mrms_hourly_correction_factor_fulres>=3.403e+37,
                                                            x = 0, y = xds_mrms_hourly_correction_factor_fulres)
    # if stageiv is non zero and mrms is zero, fill with stage iv precip intensities
    ### create stage iv data with same res as mrms
    xds_stageiv_to_mrms= spatial_resampling(ds_stageiv, ds_mrms_hourly, "latitude", "longitude")
    # replace fill values with 0 (same as above)
    xds_stageiv_to_mrms = xr.where(xds_stageiv_to_mrms>=3.403e+37, x = 0, y = xds_stageiv_to_mrms)
    ### create dataset with stageiv rainfall intensities at the indieces where condition is true and 0 everywhere else
    xds_stage_iv_where_mrms_is_0_and_stageiv_is_not = xr.where((xds_stageiv_to_mrms != 0) & (ds_mrms_hourly == 0),
                                                            x = xds_stageiv_to_mrms, y = 0)
    ### upsample to the full MRMS resolution
    #### assuming the dataset is in PROCEDING time interval, forward fill missing values
    target_index = pd.to_datetime(ds_mrms.time.values)
    xds_stage_iv_where_mrms_is_0_and_stageiv_is_not = xds_stage_iv_where_mrms_is_0_and_stageiv_is_not.reindex(dict(time = target_index)).ffill(dim="time")
    ### upsample bias correction to full res data
    xds_correction_to_mrms = xds_mrms_hourly_correction_factor_fulres.reindex(dict(time = target_index)).ffill(dim="time")
    ### apply correction factor
    xds_mrms_biascorrected = ds_mrms * xds_correction_to_mrms
    ### fill in with stageIV data where mrms data is missing
    xds_mrms_biascorrected_filled = xds_mrms_biascorrected +  xds_stage_iv_where_mrms_is_0_and_stageiv_is_not
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
    tot_rain_mrms_corrected = (xds_mrms_biascorrected_filled.mean("time")*24).rainrate.sum().values
    tot_rain_mrms_uncorrected = (ds_mrms.mean("time")*24).rainrate.sum().values
    tot_rain_stageiv = (xds_stageiv_to_mrms.mean("time")*24).rainrate.sum().values
    print("Fraction of domain-wide rainfall totals:")
    print("Bias corrected MRMS data over stage iv data: {}".format(tot_rain_mrms_corrected/tot_rain_stageiv))
    print("UNbias corrected MRMS data over stage iv data: {}".format(tot_rain_mrms_uncorrected/tot_rain_stageiv))
    #### comparing using stageiv resolution
    ### resample bias corrected mrms data to stage iv resolution to compare
    # xds_mrms_biascorrected_to_stageiv= spatial_resampling(xds_mrms_biascorrected_filled, ds_stageiv, "latitude", "longitude")
    # tot_rain_mrms_corrected = (xds_mrms_biascorrected_to_stageiv.mean("time")*24).rainrate.sum().values
    # tot_rain_mrms_uncorrected = (xds_mrms_biascorrected_to_stageiv.mean("time")*24).rainrate_uncorrected.sum().values
    # tot_rain_stageiv = (ds_stageiv.mean("time")*24).rainrate.sum().values
    # print("Fraction of domain-wide rainfall totals:")
    # print("Bias corrected MRMS data over stage iv data: {}".format(tot_rain_mrms_corrected/tot_rain_stageiv))
    # print("UNbias corrected MRMS data over stage iv data: {}".format(tot_rain_mrms_uncorrected/tot_rain_stageiv))
    xds_mrms_biascorrected_filled.attrs["domainwide_totals_CORRECTED_mrms_over_stageiv"] = tot_rain_mrms_corrected/tot_rain_stageiv
    xds_mrms_biascorrected_filled.attrs["domainwide_totals_uncorrected_mrms_over_stageiv"] = tot_rain_mrms_uncorrected/tot_rain_stageiv
    xds_mrms_biascorrected_filled.attrs["correction_factor_quantile_cutoff"] = target_quant
    xds_mrms_biascorrected_filled.attrs["correction_factor_cutoff"] = quant_correction
    return xds_mrms_biascorrected_filled, xds_correction_to_mrms, xds_stage_iv_where_mrms_is_0_and_stageiv_is_not

# xds_mrms_biascorrected_filled= bias_correct_and_fill_mrms(ds_mrms, ds_stageiv)
#%%
try:
    ds_mrms = xr.open_dataset(fl_in_nc)
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
    # print("bm 1")
    # select subset based on the extents of the transposition domain
    # the +360 is to convert from degrees west to degrees east; the + or - 0.05 is to buffer the selction by 5 gridcells assuming 0.01 degree grid
    ds_mrms = clip_ds_to_transposition_domain(ds_mrms, gdf_transdomain)
    columns.append("stageiv_available_for_bias_correction")
    if stageiv_data_available_for_bias_correction:
        values.append(True)
        ds_stageiv = xr.open_dataset(f_nc_stageiv)
        # clean up stageiv
        ds_stageiv['outlat'] = ds_stageiv.latitude.values
        ds_stageiv['outlon'] = ds_stageiv.longitude.values+360
        ds_stageiv = ds_stageiv.drop_vars("latitude")
        ds_stageiv = ds_stageiv.drop_vars("longitude")
        ds_stageiv = ds_stageiv.drop_vars("infilled")
        ds_stageiv = ds_stageiv.rename({"outlat":"latitude", "outlon":"longitude"})
        ds_stageiv = clip_ds_to_transposition_domain(ds_stageiv, gdf_transdomain)
        ds_mrms_biascorrected_filled,__,__ = bias_correct_and_fill_mrms(ds_mrms, ds_stageiv)
        columns.append("domainwide_totals_CORRECTED_mrms_over_stageiv")
        values.append(ds_mrms_biascorrected_filled.attrs["domainwide_totals_CORRECTED_mrms_over_stageiv"])
        columns.append("domainwide_totals_uncorrected_mrms_over_stageiv")
        values.append(ds_mrms_biascorrected_filled.attrs["domainwide_totals_uncorrected_mrms_over_stageiv"])
        columns.append("correction_factor_quantile_cutoff")
        values.append(ds_mrms_biascorrected_filled.attrs["correction_factor_quantile_cutoff"])
        columns.append("correction_factor_cutoff")
        values.append(ds_mrms_biascorrected_filled.attrs["correction_factor_cutoff"])
    else:
        values.append(False)
    # print("bm 4")
    df_input_dataset_attributes = pd.DataFrame([values], columns=columns)
    # print("bm 5")
    # df_input_dataset_attributes = pd.DataFrame(ds.attrs, index = [0]) # the attributes are the columns, index is the filepath to the netcdf
    # df_input_dataset_attributes['filepath'] = [ds_mrms.encoding['source']]
    # print("bm 6")
except Exception as e:
    print("The following error was encountered:")
    print(e)
    performance["loading_netcdf_errors"]  = e
    performance["problem_loading_netcdf"] = True

# tstep = ds.attrs["time_step"]
if performance["problem_loading_netcdf"] == False:
    # verify the full day has coverage
    tstep_min = pd.to_timedelta(ds_mrms.attrs["time_step"]).total_seconds() / 60
    num_tsteps = ds_mrms_biascorrected_filled.coords["time"].shape[0]
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
        t_idx_1min = pd.date_range(ds_mrms_biascorrected_filled.time.values[0], periods = 24*60, freq='1min')
        ds_1min = ds_mrms_biascorrected_filled.reindex(dict(time = t_idx_1min)).ffill(dim="time")
        da_target = ds_1min.resample(time = "{}Min".format(target_tstep)).mean()
        # performance["problems_resampling"] = False
    else:
        da_target = ds_mrms_biascorrected_filled
    performance["problem_exporting_zarr"] = False
    performance["to_zarr_errors"] = "None"
    try:
        da_target.to_zarr(fl_out_zarr, mode="w")
        ds_from_zarr = xr.open_zarr(store=fl_out_zarr)
    except Exception as e:
        performance["to_zarr_errors"]  = e
        performance["problem_exporting_zarr"] = True
    # export to netcdf
    performance["problem_exporting_netcdf"] = False
    performance["to_netcdf_errors"] = "None"
    try:
        # ds_from_zarr.to_netcdf(fl_out_nc, encoding= {"rainrate":{"zlib":True},"rainrate_uncorrected":{"zlib":True},
        #                        "mrms_bias_correction_multiplier":{"zlib":True},"stageiv_fillvals_where_mrms_is_0_and_stageiv_is_not":{"zlib":True}})
        ds_from_zarr.to_netcdf(fl_out_nc, encoding= {"rainrate":{"zlib":True}})
        shutil.rmtree(fl_out_zarr)
    except Exception as e:
        print("Exporting netcdf dataset failed due to error: {}".format(e))
        performance["to_netcdf_errors"]  = e
        performance["problem_exporting_netcdf"] = True

# export performance dictionary to a csv
time_elapsed_min = round((time.time() - start_time) / 60, 2)
performance["time_elapsed_min"] = time_elapsed_min
df = pd.DataFrame(performance, index = [1])
df.to_csv(fl_out_csv)
if performance["problem_loading_netcdf"] == False:
    df_input_dataset_attributes.to_csv(fl_out_csv_qaqc)
print("script finished")
# %%
