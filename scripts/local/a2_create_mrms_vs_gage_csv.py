#%% Import libraries and set filepaths
import sys
import pandas as pd
import pytz
import glob
import dask
dask.config.set(**{'array.slicing.split_large_chunks': False})

inches_per_mm = 1/25.4
f_in_csvs = "D:/mrms_processing/data/mrms_for_rainyday_subset_norfolk_csvs/*.csv"
fl_events = "D:/Dropbox/_GradSchool/_norfolk/norfolk_ffa/output_csvs_from_code/rainfall_events_all.csv"
f_out_csv = "D:/mrms_processing/out_csvs/a_mrms_and_hrsd_event_data.csv"
fl_csv_prior_processing = "D:/mrms_processing/data/data/mrms_at_gages_prior_processing_method/2015-1-1_to_2021-9-22_mrms_at_gages.csv"
#%% load mrms data from csvs
all_files = glob.glob(f_in_csvs)
# load data of mrms gridcells overlapping HRSD gages
df_mrms_long = pd.concat((pd.read_csv(f, index_col=0, parse_dates=["time"]) for f in all_files), ignore_index=True)

# assign timezone to time
dt = df_mrms_long.time
dti = pd.DatetimeIndex(dt, tz="utc")
df_mrms_long.time = dti
#%% load HRSD data
gage_ids = df_mrms_long.overlapping_gage_id.unique()

df_events_wide = pd.read_csv(fl_events, parse_dates=["Time"])
dti = pd.DatetimeIndex(df_events_wide.Time)
if dti.tz != pytz.utc:
    sys.exit("The timezone of the gage data is either non-UTC or unassigned. Exiting script...")

df_events_long = pd.melt(df_events_wide, id_vars = ["Time", "event_id"], value_vars = gage_ids,
                        var_name = "gage_id", value_name = "gage_rainfall")
#%% convert mrms data to a proceeding 15-minute timestamp (note that it is already a proceeding time interval)
df_mrms_long = df_mrms_long.sort_values(["overlapping_gage_id", "time"])

## convert to 1 minute 
df_mrms_long_1min = df_mrms_long.set_index("time")
df_mrms_long_1min = df_mrms_long_1min.groupby("overlapping_gage_id")
df_mrms_long_1min = df_mrms_long_1min.resample("1min").bfill()
df_mrms_long_1min = df_mrms_long_1min.droplevel("overlapping_gage_id")

# convert to 15 minute
df_mrms_long_15min = df_mrms_long_1min.groupby(["overlapping_gage_id", "mrms_lat", "mrms_long"])
df_mrms_long_15min = df_mrms_long_15min.resample("15min").mean()
df_mrms_long_15min = df_mrms_long_15min.drop(labels=["mrms_lat", "mrms_long"], axis="columns")

# convert from mm/hr to inches
df_mrms_long_15min["mrms_precip_in"] = df_mrms_long_15min.loc[:, "precip_mm_per_hour"] * 15/60 * inches_per_mm # mm/hr * 15/60 hrs * in/mm
#%% merge the dataframe 
df_out = pd.merge(df_mrms_long_15min, df_events_long, how = 'left', left_on = ['time', 'overlapping_gage_id'], right_on = ['Time', 'gage_id'])

# drop timesteps where HRSD events and MRMS timesteps do not overlap
df_out.dropna(subset="event_id", inplace=True)

df_out =  df_out.reset_index(drop=True)

df_out = df_out.rename(columns={"Time":"time", "gage_rainfall":"gage_precip_in"})

df_out = df_out.loc[:, ['time', 'gage_id', 'event_id', 'mrms_precip_in', 'gage_precip_in']]

df_out = df_out.astype({"event_id":int})

#%% export csv
df_out.to_csv(f_out_csv)
print("finished exporting csv")