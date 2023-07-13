#%% Import libraries
import shutil
import pandas as pd
import sys
from glob import glob

target_tstep = 5

chnk_sz = "5000MB"

performance = {}
#%% work
# fldr_in_nc_day = "D:/Dropbox/_GradSchool/_norfolk/highres-radar-rainfall-processing/data/mrms_nc_preciprate_fullres_dailyfiles/"
# fldr_out_nc = "D:/Dropbox/_GradSchool/_norfolk/highres-radar-rainfall-processing/data/_scratch/".format(target_tstep)
# in_date = "20190306"

# fldr_in_nc_day = "/project/quinnlab/dcl3nd/norfolk/highres-radar-rainfall-processing/data/mrms_nc_preciprate_fullres_dailyfiles/"
fldr_out_nc = "/scratch/dcl3nd/highres-radar-rainfall-processing/out_fullres_dailyfiles_consolidated/"
# fldr_out_zarr = "/project/quinnlab/dcl3nd/norfolk/highres-radar-rainfall-processing/data/_scratch/zarrs/"
fldr_out_csv = "/project/quinnlab/dcl3nd/norfolk/highres-radar-rainfall-processing/data/_scratch/csv/"
#%% end work

# folders (with proceeding fwd slash)
# fldr_in_nc_day = str(sys.argv[1]) # ${assar_dirs[out_fullres_dailyfiles]} # "/project/quinnlab/dcl3nd/norfolk/highres-radar-rainfall-processing/data/mrms_nc_preciprate_fullres_dailyfiles/"
fldr_out_nc = str(sys.argv[2]) # ${assar_dirs[out_fullres_dailyfiles_consolidated]} # "/scratch/dcl3nd/highres-radar-rainfall-processing/out_fullres_dailyfiles_consolidated/"
# fldr_out_zarr = str(sys.argv[3]) # ${assar_dirs[scratch_zarrs]} # "/project/quinnlab/dcl3nd/norfolk/highres-radar-rainfall-processing/data/_scratch/zarrs/"
fldr_in_csv = str(sys.argv[4]) # ${assar_dirs[scratch_zarrs]} # "/project/quinnlab/dcl3nd/norfolk/highres-radar-rainfall-processing/data/_scratch/csv/"

#%% netcdf 
# fl_in_nc = fldr_in_nc_day +"{}.nc".format(in_date)
# fl_out_nc = fldr_out_nc +"{}.nc".format(in_date)
# fl_out_zarr = fldr_out_zarr +"{}.zarr".format(in_date)
fl_out_csv = fldr_out_csv +"da2_resampling_{}.csv".format("*") # must match pattern in script da2
lst_f_csvs = glob(fl_out_csv)

lst_dfs = []
for f in lst_f_csvs:
    lst_dfs.append(pd.read_csv(f, index_col = 0))

df = pd.concat(lst_dfs, ignore_index = True)

df.to_csv(fldr_out_nc+"_da3_resampling_performance.csv")

#%% work
# df[(df.problem_loading_netcdf == False) & (df.current_tstep_different_than_target == True)].drop(["problem_loading_netcdf", "loading_netcdf_errors", "duration_h", "problem_with_duration", "current_tstep_different_than_target"], axis = 1)
