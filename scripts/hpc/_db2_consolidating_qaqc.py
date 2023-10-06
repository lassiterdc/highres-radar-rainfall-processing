#%%
import shutil
import pandas as pd
import sys
from glob import glob
#%% # inputs

f_out_csv = str(sys.argv[1]) + "db_consolidating_tseps_{}.csv".format("*")
fldr_out_nc = str(sys.argv[2]) + "/_qaqc_mrms_nc_preciprate_hourly_and_daily_dailyfiles.csv"
lst_f_csvs = glob(f_out_csv)

lst_dfs = []
for f in lst_f_csvs:
    lst_dfs.append(pd.read_csv(f))

df = pd.concat(lst_dfs, ignore_index = True)

df.to_csv(fldr_out_nc+"_da3_resampling_performance.csv")