#%% Import libraries
import shutil
import pandas as pd
import sys
from glob import glob
import xarray as xr
import time

start_time = time.time()

target_tstep = 5

chnk_sz = "5000MB"

performance = {}
#%% work
fldr_out_nc = "/scratch/dcl3nd/highres-radar-rainfall-processing/data/mrms_nc_preciprate_fullres_dailyfiles_constant_tstep/"
fldr_csvs = "/scratch/dcl3nd/highres-radar-rainfall-processing/data/_scratch/csv/"
#%% end work
fldr_out_nc = str(sys.argv[1]) # ${assar_dirs[out_fullres_dailyfiles_consolidated]} # "/scratch/dcl3nd/highres-radar-rainfall-processing/out_fullres_dailyfiles_consolidated/"
fldr_csvs = str(sys.argv[2])

fl_da2_csv = fldr_csvs +"da2_resampling_{}.csv".format("*") # must match pattern in script da2

# qaqc of resampling
lst_f_csvs = glob(fl_da2_csv)

lst_dfs = []
for f in lst_f_csvs:
    lst_dfs.append(pd.read_csv(f, index_col = 0))

df = pd.concat(lst_dfs, ignore_index = True)

df_successes = df[df.problem_exporting_netcdf == False]

years = pd.to_datetime(df_successes.date,format="%Y%m%d").dt.year.unique()

df.to_csv(fldr_out_nc+"_da3_resampling_performance.csv")

# qaqc of dataset as a whole
fl_csv_qaqc = fldr_csvs +"qaqc_of_daily_fullres_data_{}.csv".format("*") # must match pattern in script da2
lst_f_csvs = glob(fl_csv_qaqc)
lst_dfs = []
for f in lst_f_csvs:
    lst_dfs.append(pd.read_csv(f, index_col = 0))
df = pd.concat(lst_dfs, ignore_index = False)
df.to_csv(fldr_out_nc+"_da3_qaqc_fullres_nc_dataset.csv")

# consolidating qaqc statistics in the netcdfs
for year in years:
    lst_f_netcdfs = glob(fldr_out_nc + "{}*.nc".format(year)) # 2 is there to ensure that the processed data are the only netcdfs that are included
    # DCL WORK
    # lst_f_netcdfs = lst_f_netcdfs[0:366]
    # END DCL WORK
    lst_f_netcdfs.sort()
    lst_ds_qaqc = []
    for f_nc in lst_f_netcdfs:
        ds_day = xr.open_dataset(f_nc)
        lst_das = []
        date = ds_day.time.values[0]
        # print(date)
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
    fl_nc_qaqc_out = fldr_out_nc + "_qaqc_of_resampled_data_{}.nc".format(year)
    ds_qaqc.to_netcdf(fl_nc_qaqc_out, encoding=d_encoding, engine="netcdf4")
    time_elapsed_min = round((time.time() - start_time) / 60, 2)
    print("Exported qaqc file for year {}. Total elapsed time (min): {}. File: ".format(year, time_elapsed_min, fl_nc_qaqc_out))

print("Finishe running script in {} minutes".format(time_elapsed_min))



    