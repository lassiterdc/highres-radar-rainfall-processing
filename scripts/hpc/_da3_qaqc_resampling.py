#%% Import libraries
import shutil
import pandas as pd
import sys
from glob import glob
import xarray as xr
import time
from __utils import *
import os

start_time = time.time()

#%% work
fldr_out_nc = "/scratch/dcl3nd/highres-radar-rainfall-processing/data/mrms_nc_preciprate_fullres_dailyfiles_constant_tstep/"
fldr_csvs = "/scratch/dcl3nd/highres-radar-rainfall-processing/data/_scratch/csv/"
fl_nc_qaqc_out_all = fldr_out_nc + "_qaqc_of_resampled_data.nc"
fl_da2_csv = fldr_csvs +"da2_resampling_{}.csv".format("*") # must match pattern in script da2
# use_previous_qaqc_netcdf = True
#%% end work
fldr_out_nc = str(sys.argv[1]) # ${assar_dirs[out_fullres_dailyfiles_consolidated]} # "/scratch/dcl3nd/highres-radar-rainfall-processing/out_fullres_dailyfiles_consolidated/"
fldr_csvs = str(sys.argv[2])

fl_nc_qaqc_out_all = fldr_out_nc + "_qaqc_of_resampled_data.nc"
fl_nc_qaqc_out_all_st4_res = fl_nc_qaqc_out_all.split(".nc")[0] + "_st4_res.nc"
fl_da2_csv = fldr_csvs +"da2_resampling_{}.csv".format("*") # must match pattern in script da2

# qaqc of resampling
lst_f_csvs = glob(fl_da2_csv)

lst_dfs = []
for f in lst_f_csvs:
    lst_dfs.append(pd.read_csv(f, index_col = 0))

df = pd.concat(lst_dfs, ignore_index = True)

df_successes = df[df.problem_exporting_netcdf == False]

df.to_csv(fldr_out_nc+"_da3_resampling_performance.csv")

# qaqc of dataset as a whole
fl_csv_qaqc = fldr_csvs +"qaqc_of_daily_fullres_data_{}.csv".format("*") # must match pattern in script da2
lst_f_csvs = glob(fl_csv_qaqc)
lst_dfs = []
for f in lst_f_csvs:
    lst_dfs.append(pd.read_csv(f, index_col = 0))
df = pd.concat(lst_dfs, ignore_index = False)
df.to_csv(fldr_out_nc+"_da3_qaqc_fullres_nc_dataset.csv")

# consolidating qaqc statistics into monthly netcdf files
df_dates = pd.to_datetime(df_successes.date,format="%Y%m%d").sort_values().astype(str).str.split("-", expand = True)
df_dates.columns = ["year", "month", "day"]
df_yearmonths = df_dates.loc[:, ["year", "month"]].drop_duplicates()
df_yearmonths = df_yearmonths.reset_index(drop=True)

for id, row in df_yearmonths.iterrows():
    bm_time = time.time()
    year = row.year
    month = row.month
    lst_f_netcdfs = glob(fldr_out_nc + "{}{}*.nc".format(year, month))
    lst_f_netcdfs.sort()
    lst_ds_qaqc = []
    for f_nc in lst_f_netcdfs:
        ds_day = xr.open_dataset(f_nc)
        lst_das = []
        date = ds_day.time.values[0]
        for dvar in ds_day.data_vars:
            include = True
            for coord in ds_day[dvar].coords: # if time is one of the coordinates, do not include it
                if "time" in coord:
                    include = False
            if include:
                da_to_include = ds_day[dvar].reset_coords()[dvar] # isolate just the dataarray and the used coordinates
                # add as a coordinate the date
                da_to_include = da_to_include.assign_coords({"date":date})
                lst_das.append(da_to_include)
        ds_day_qaqc = xr.Dataset({da.name: da for da in lst_das})
        lst_ds_qaqc.append(ds_day_qaqc)
    ds_qaqc = xr.combine_nested(lst_ds_qaqc, concat_dim = "date")
    d_encoding = {}
    for da_name in ds_qaqc.data_vars:
        d_encoding[da_name] = {"zlib":True}
    fl_scratch_nc_qaqc = fldr_scratch_zarr + "da3_qaqc_{}{}.nc".format(year, month)
    ds_qaqc.load().to_netcdf(fl_scratch_nc_qaqc, encoding=d_encoding, engine="h5netcdf")
    time_elapsed_min = round((time.time() - bm_time) / 60, 2)
    print("Exported temporary netcdf qaqc file for year and month {}-{}. Time to export: {}.".format(year, month, time_elapsed_min))

# load one year's worth of zarr files and consolidate then
    
# also create version with stage iv dimensions
ds_qaqc_all = xr.open_dataset(fl_nc_qaqc_out_all, engine="h5netcdf", chunks = dict(date = 22))
# load stage iv data
lst_f_st4 = glob(fldr_nc_stageiv + "*/*nc")
lst_f_st4.sort()
f_latest_st4 = lst_f_st4[-1]
ds_stageiv = xr.open_dataset(f_latest_st4)
ds_stageiv = process_dans_stageiv(ds_stageiv)
ds_stageiv_subset = clip_ds_to_another_ds(ds_stageiv, ds_qaqc_all)
lst_ncs_year_st4_res = []

lst_ncs_year = []
for year in df_yearmonths.year.unique():
    bm_time = time.time()
    lst_f_scratch_nc = glob(fldr_scratch_zarr + "da3_qaqc_{}*.nc".format(year))
    ds_qaqc_year = xr.open_mfdataset(lst_f_scratch_nc, engine = "h5netcdf")
    fl_nc_qaqc_out = fldr_out_nc + "_qaqc_of_resampled_data_{}.nc".format(year)
    d_encoding = {}
    for da_name in ds_qaqc_year.data_vars:
        d_encoding[da_name] = {"zlib":True}
    ds_qaqc_year.to_netcdf(fl_nc_qaqc_out, encoding=d_encoding, engine="h5netcdf")
    time_elapsed_min = round((time.time() - bm_time) / 60, 2)
    lst_ncs_year.append(fl_nc_qaqc_out)
    print("Exported netcdf qaqc file for year {}. Time to export: {}.".format(year, time_elapsed_min))
    bm_time2 = time.time()
    ds_qaqc_year_to_st4 = spatial_resampling(ds_qaqc_year, ds_stageiv_subset, "latitude", "longitude", missingfillval = 0)
    f_out_resampled = fl_nc_qaqc_out.split(".nc")[0] + "_st4_res.nc"
    for da_name in ds_qaqc_year_to_st4.data_vars:
        d_encoding[da_name] = {"zlib":True}
    ds_qaqc_year_to_st4.to_netcdf(f_out_resampled, encoding=d_encoding, engine="h5netcdf")
    lst_ncs_year_st4_res.append(f_out_resampled)
    time_elapsed_min = round((time.time() - bm_time2) / 60, 2)
    print("Additional time to export netcdf consolidated to stage iv resolution: {}.".format(time_elapsed_min))
    

# finally combine them all into a single netcdf
## mrms spatial resolution
bm_time = time.time()
ds_qaqc_all = xr.open_mfdataset(lst_ncs_year, engine = "h5netcdf")
d_encoding = {}
for da_name in ds_qaqc_all.data_vars:
    d_encoding[da_name] = {"zlib":True}
ds_qaqc_all.to_netcdf(fl_nc_qaqc_out_all, encoding=d_encoding, engine="h5netcdf")
time_elapsed_min = round((time.time() - bm_time) / 60, 2)
print("Exported netcdf qaqc file for entire dataset. Time to export: {}.".format(time_elapsed_min))
print("filepath exported: {}".format(fl_nc_qaqc_out_all))
# scratch files once netcdf has been created
for f in (lst_f_scratch_nc+lst_ncs_year):
    # shutil.rmtree(f)
    os.remove(f)   

## stage iv spatial resolution
bm_time = time.time()
ds_qaqc_all = xr.open_mfdataset(lst_ncs_year_st4_res, engine = "h5netcdf")
d_encoding = {}
for da_name in ds_qaqc_all.data_vars:
    d_encoding[da_name] = {"zlib":True}

ds_qaqc_all.to_netcdf(fl_nc_qaqc_out_all_st4_res, encoding=d_encoding, engine="h5netcdf")
time_elapsed_min = round((time.time() - bm_time) / 60, 2)
print("Exported netcdf qaqc file for entire dataset. Time to export: {}.".format(time_elapsed_min))
print("filepath exported: {}".format(fl_nc_qaqc_out_all_st4_res))
# scratch files once netcdf has been created
for f in (lst_ncs_year_st4_res):
    # shutil.rmtree(f)
    os.remove(f)  

