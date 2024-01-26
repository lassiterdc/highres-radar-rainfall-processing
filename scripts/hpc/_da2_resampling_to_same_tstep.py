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

target_tstep = 5

chnk_sz = "1000MB"

performance = {}
#%% work
# fldr_nc_fullres_daily = "D:/Dropbox/_GradSchool/_norfolk/highres-radar-rainfall-processing/data/mrms_nc_preciprate_fullres_dailyfiles/"
# fldr_nc_fullres_daily_constant_tstep = "D:/Dropbox/_GradSchool/_norfolk/highres-radar-rainfall-processing/data/_scratch/".format(target_tstep)
# in_date = "20190306"
# in_date = "20160410"





#%% end work

# folders (with proceeding fwd slash)
in_date = str(sys.argv[1]) # YYYYMMDD
fldr_nc_fullres_daily = str(sys.argv[2]) # ${assar_dirs[out_fullres_dailyfiles]} # "/project/quinnlab/dcl3nd/norfolk/highres-radar-rainfall-processing/data/mrms_nc_preciprate_fullres_dailyfiles/"
fldr_nc_fullres_daily_constant_tstep = str(sys.argv[3]) # ${assar_dirs[out_fullres_dailyfiles_consolidated]} # "/scratch/dcl3nd/highres-radar-rainfall-processing/out_fullres_dailyfiles_consolidated/"
fldr_scratch_zarr = str(sys.argv[4]) # ${assar_dirs[scratch_zarrs]} # "/project/quinnlab/dcl3nd/norfolk/highres-radar-rainfall-processing/data/_scratch/zarrs/"
fldr_scratch_csv = str(sys.argv[5]) # ${assar_dirs[scratch_zarrs]} # "/project/quinnlab/dcl3nd/norfolk/highres-radar-rainfall-processing/data/_scratch/csv/"
f_shp_sst_transom = str(sys.argv[6]) # ${assar_dirs[shp_transposition_domain]} # "/project/quinnlab/dcl3nd/norfolk/stormy/stochastic_storm_transposition/norfolk/transposition_domain/norfolk_trans_dom_4326.shp"

performance["date"] = in_date

f_out_export_perf = fldr_scratch_zarr + "_export_stats_{}.csv".format(in_date)
#%% netcdf 
fl_in_nc = fldr_nc_fullres_daily +"{}.nc".format(in_date)
fl_out_nc = fldr_nc_fullres_daily_constant_tstep +"{}.nc".format(in_date)
fl_out_zarr = fldr_scratch_zarr +"{}.zarr".format(in_date)
fl_out_csv = fldr_scratch_csv +"da2_resampling_{}.csv".format(in_date)

performance["problem_loading_netcdf"] = False
performance["loading_netcdf_errors"]  = "None"
try:
    ds = xr.open_dataset(fl_in_nc)
    # select subset based on the extents of the transposition domain
    gdf_transdomain = gp.read_file(f_shp_sst_transom)
    transdom_bounds = gdf_transdomain.bounds
    # the +360 is to convert from degrees west to degrees east; the + or - 0.05 is to buffer the selction by 5 gridcells assuming 0.01 degree grid
    ds = ds.where((ds.latitude >= float(transdom_bounds.miny-.05)) & (ds.latitude <= float(transdom_bounds.maxy+.05)) & (ds.longitude >= float(transdom_bounds.minx+360-.05)) & (ds.longitude <= float(transdom_bounds.maxx+360+.05)), drop = True)
except Exception as e:
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
        performance["problems_resampling"] = True
        t_idx_1min = pd.date_range(ds.time.values[0], periods = 24*60, freq='1min')
        ds_1min = ds.reindex(dict(time = t_idx_1min)).ffill(dim="time")
        da_target = ds_1min.resample(time = "{}Min".format(target_tstep)).mean()
        performance["problems_resampling"] = False
        # WORK
        # import numpy as np
        # da_target = da_target.isel(dict(time = np.arange(100, 110)))
        # END WORK
        # export to zarr
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
    # else: # if the timestep is already correct, just copy the file
    #     shutil.copyfile(fl_in_nc, fl_out_nc)

# export performance dictionary to a csv
time_elapsed_min = round((time.time() - start_time) / 60, 2)
performance["time_elapsed_min"] = time_elapsed_min
df = pd.DataFrame(performance, index = [1])
df.to_csv(fl_out_csv)