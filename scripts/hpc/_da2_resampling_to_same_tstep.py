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

target_tstep = 5

chnk_sz = "1000MB"

performance = {}
#%% work

# in_date = "20210110"
# fldr_nc_fullres_daily = "/scratch/dcl3nd/highres-radar-rainfall-processing/data/mrms_nc_preciprate_fullres_dailyfiles/"
# fldr_nc_fullres_daily_constant_tstep = "/scratch/dcl3nd/highres-radar-rainfall-processing/out_fullres_dailyfiles_consolidated/"
# fldr_scratch_zarr = "/scratch/dcl3nd/highres-radar-rainfall-processing/data/_scratch/zarrs/"
# fldr_scratch_csv = "/scratch/dcl3nd/highres-radar-rainfall-processing/data/_scratch/csv/"
# f_shp_sst_transom = "/scratch/dcl3nd/stormy/stochastic_storm_transposition/norfolk/transposition_domain/norfolk_trans_dom_4326.shp"


#%% end work



# folders (with proceeding fwd slash)
in_date = str(sys.argv[1]) # YYYYMMDD
fldr_nc_fullres_daily = str(sys.argv[2]) # ${assar_dirs[out_fullres_dailyfiles]} # "/project/quinnlab/dcl3nd/norfolk/highres-radar-rainfall-processing/data/mrms_nc_preciprate_fullres_dailyfiles/"
fldr_nc_fullres_daily_constant_tstep = str(sys.argv[3]) # ${assar_dirs[out_fullres_dailyfiles_consolidated]} # "/scratch/dcl3nd/highres-radar-rainfall-processing/out_fullres_dailyfiles_consolidated/"
fldr_scratch_zarr = str(sys.argv[4]) # ${assar_dirs[scratch_zarrs]} # "/project/quinnlab/dcl3nd/norfolk/highres-radar-rainfall-processing/data/_scratch/zarrs/"
fldr_scratch_csv = str(sys.argv[5]) # ${assar_dirs[scratch_zarrs]} # "/project/quinnlab/dcl3nd/norfolk/highres-radar-rainfall-processing/data/_scratch/csv/"
f_shp_sst_transom = str(sys.argv[6]) # ${assar_dirs[shp_transposition_domain]} # "/project/quinnlab/dcl3nd/norfolk/stormy/stochastic_storm_transposition/norfolk/transposition_domain/norfolk_trans_dom_4326.shp"

performance["date"] = in_date

# f_out_export_perf = fldr_scratch_zarr + "_export_stats_{}.csv".format(in_date)
#%% netcdf 
fl_in_nc = fldr_nc_fullres_daily +"{}.nc".format(in_date)
fl_out_nc = fldr_nc_fullres_daily_constant_tstep +"{}.nc".format(in_date)
fl_out_zarr = fldr_scratch_zarr +"{}.zarr".format(in_date)
fl_out_csv = fldr_scratch_csv +"da2_resampling_{}.csv".format(in_date)
fl_out_csv_qaqc = fldr_scratch_csv +"qaqc_of_daily_fullres_data_{}.csv".format(in_date)

performance["problem_loading_netcdf"] = False
performance["loading_netcdf_errors"]  = "None"
try:
    ds = xr.open_dataset(fl_in_nc)
    # create a single row dataset with netcdf attributes
    columns = []
    values = []
    for key in ds.attrs:
        if isinstance(ds.attrs[key], dict):
            for key2 in ds.attrs[key]:
                columns.append(key2)
                values.append(str(ds.attrs[key][key2]))
        if key == 'warnings':
            input_string = ds.attrs[key]
            result_dict = ast.literal_eval(input_string)
            for key2 in result_dict:
                columns.append(key2)
                values.append(str(result_dict[key2]))
        else:
            columns.append(key)
            values.append(str(ds.attrs[key]))
    df_input_dataset_attributes = pd.DataFrame([values], columns=columns)
    # df_input_dataset_attributes = pd.DataFrame(ds.attrs, index = [0]) # the attributes are the columns, index is the filepath to the netcdf
    df_input_dataset_attributes['filepath'] = [ds.encoding['source']]
    # select subset based on the extents of the transposition domain
    gdf_transdomain = gp.read_file(f_shp_sst_transom)
    transdom_bounds = gdf_transdomain.bounds
    # the +360 is to convert from degrees west to degrees east; the + or - 0.05 is to buffer the selction by 5 gridcells assuming 0.01 degree grid
    ds = ds.where((ds.latitude >= float(transdom_bounds.miny.iloc[0]-.05)) & (ds.latitude <= float(transdom_bounds.maxy.iloc[0]+.05)) & (ds.longitude >= float(transdom_bounds.minx.iloc[0]+360-.05)) & (ds.longitude <= float(transdom_bounds.maxx.iloc[0]+360+.05)), drop = True)
except Exception as e:
    print("The following error was encountered:")
    print(e)
    performance["loading_netcdf_errors"]  = e
    performance["problem_loading_netcdf"] = True

# tstep = ds.attrs["time_step"]
if performance["problem_loading_netcdf"] == False:
    # verify the full day has coverage
    tstep_min = pd.to_timedelta(ds.attrs["time_step"]).total_seconds() / 60
    num_tsteps = ds.coords["time"].shape[0]
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
        t_idx_1min = pd.date_range(ds.time.values[0], periods = 24*60, freq='1min')
        ds_1min = ds.reindex(dict(time = t_idx_1min)).ffill(dim="time")
        da_target = ds_1min.resample(time = "{}Min".format(target_tstep)).mean()
        # performance["problems_resampling"] = False
    else:
        da_target = ds
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
        ds_from_zarr.to_netcdf(fl_out_nc, encoding= {"rainrate":{"zlib":True}})
        shutil.rmtree(fl_out_zarr)
    except Exception as e:
        # print("Simulation failed due to error: {}".format(e))
        performance["to_netcdf_errors"]  = e
        performance["problem_exporting_netcdf"] = True

# export performance dictionary to a csv
time_elapsed_min = round((time.time() - start_time) / 60, 2)
performance["time_elapsed_min"] = time_elapsed_min
df = pd.DataFrame(performance, index = [1])
df.to_csv(fl_out_csv)
if performance["problem_loading_netcdf"] == False:
    df_input_dataset_attributes.to_csv(fl_out_csv_qaqc)
# %%
