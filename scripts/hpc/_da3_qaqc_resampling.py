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

# fldr_in_nc_day = "/scratch/dcl3nd/highres-radar-rainfall-processing/data/mrms_nc_preciprate_fullres_dailyfiles/"
# fldr_out_nc = "/scratch/dcl3nd/highres-radar-rainfall-processing/data/mrms_nc_preciprate_fullres_dailyfiles_constant_tstep/"
# fldr_out_zarr = "/scratch/dcl3nd/highres-radar-rainfall-processing/data/_scratch/zarrs/"
# fldr_out_csv = "/scratch/dcl3nd/highres-radar-rainfall-processing/data/_scratch/csv/"
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

df.to_csv(fldr_out_nc+"_da3_resampling_performance.csv")

# qaqc of dataset as a whole
fl_csv_qaqc = fldr_csvs +"qaqc_of_daily_fullres_data_{}.csv".format("*") # must match pattern in script da2
lst_f_csvs = glob(fl_csv_qaqc)
lst_dfs = []
for f in lst_f_csvs:
    lst_dfs.append(pd.read_csv(f, index_col = 0))
df = pd.concat(lst_dfs, ignore_index = False)
df.to_csv(fldr_out_nc+"_da3_qaqc_fullres_nc_dataset.csv")