#%% Import libraries and set filepaths
from operator import index
from turtle import color
import xarray as xr
import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
import itertools
import dask
dask.config.set(**{'array.slicing.split_large_chunks': False})
import shutil
import cartopy.crs as ccrs
import geopandas as gpd
import scipy.stats as st
from __filepaths import return_b_filepaths

exclude_2012_2013_and_2014 = True

inches_per_mm = 1/25.4
# f_in_nc_24h = "D:/mrms_processing/out_netcdfs/a1_mrms_at_norfolk_24h.nc"
# f_in_nc = "D:/mrms_processing/out_netcdfs/a1_mrms_at_norfolk.nc"
# f_in_csv = "D:/mrms_processing/out_csvs/a_mrms_and_hrsd_event_data.csv"
# fldr_out_plots = "D:/mrms_processing/plots/b_visualizations_of_gage_vs_mrms_data/"
# f_shp_coast = "D:/mrms_processing/data/shapefiles/composite_shoreline_clipped.shp"
# f_shp_subcatchments = "D:/mrms_processing/data/shapefiles/subcatchments.shp"
# f_shp_gages = 'D:/mrms_processing/data/gage_hrsd/rain_gages.shp'

f_nc_24h_atgages, f_nc_atgages, f_csv_mrms_and_gage_events, fldr_out_plots, f_shp_coast, f_shp_subcatchments, f_shp_gages = return_b_filepaths()
#%% analyze event totals
# loading data
df_gage_and_mrms = pd.read_csv(f_csv_mrms_and_gage_events, parse_dates=["time"], infer_datetime_format=True, index_col=0)

mrms_1d = xr.open_dataset(f_nc_24h_atgages).load()
mrms = xr.open_dataset(f_nc_atgages, chunks={"time":"1000MB"})

def create_time_mask_to_exclude_target_years(tseries):
    dti_years = pd.DatetimeIndex(tseries).year
    exclude = dti_years.where((dti_years != 2012) & (dti_years != 2013) & (dti_years != 2014))
    mask = exclude >= 0 # this yields a False at nan locations
    return mask

if exclude_2012_2013_and_2014 == True:
    mask = create_time_mask_to_exclude_target_years(df_gage_and_mrms.time)
    df_gage_and_mrms = df_gage_and_mrms[mask]

    mask = create_time_mask_to_exclude_target_years(mrms_1d.time.values)
    times_idx = (np.where(mask == True))[0]
    mrms_1d = mrms_1d.isel(time=times_idx)

    mask = create_time_mask_to_exclude_target_years(mrms.time.values)
    times_idx = (np.where(mask == True))[0]
    mrms = mrms.isel(time=times_idx)


#%% plotting using geopandas
figwidth, figheight = (7.5, 6)
# gds_states = gpd.read_file(f_shp_states)
gds_subs = gpd.read_file(f_shp_subcatchments)
gds_gages = gpd.read_file(f_shp_gages)
gds_coast = gpd.read_file(f_shp_coast)

proj = ccrs.PlateCarree()
gds_subs_trns = gds_subs.to_crs(proj)
# merge subcatchments into one
# polygon = gds_subs_trns.geometry.unary_union
# gds_subs_trns = gpd.GeoDataFrame(geometry=[polygon], crs=gds_subs_trns.crs)

# gds_states_trns = gds_states.to_crs(proj)
gds_gages_trns = gds_gages.to_crs(proj)
gds_coast_trns = gds_coast.to_crs(proj)

def convert_mm_per_hour_to_in_per_event(da):
    tseries = da.time.values
    hrs_in_event = (max(tseries) - min(tseries)) / np.timedelta64(1, 'h')
    da = da * hrs_in_event * inches_per_mm # mm/hr * hrs in event * inches_per_mm
    return da

def plot_day_of_mrms_data(gage_totals, tseries, cbar_high, title, fname, gds_subs_trns=gds_subs_trns, gds_gages_trns=gds_gages_trns, gds_coast_trns=gds_coast_trns):
    fig, ax = plt.subplots(figsize = [figwidth, figheight], subplot_kw=dict(projection=proj), dpi=300)
    # cbar_high = mrms_1d.rainrate.max() * inches_per_mm

    # modifying gages geodatabase
    gage_totals_renamed = gage_totals.reset_index().rename(columns={"gage_id":"MONITORING"})
    gds_gages_trns_with_rainfall = gds_gages_trns.merge(gage_totals_renamed, on="MONITORING")
    
    # subsetting based on time
    tseries_min = min(tseries.values)
    tseries_max = max(tseries.values)
    v_time = mrms.rainrate.time.values
    v_time = v_time[(v_time>=tseries_min) & (v_time<=tseries_max)]

    # troubleshooting 2015 not showing up 
    dti = pd.DatetimeIndex(mrms.rainrate.time.values)
    if sum(dti.year == 2015) == 0:
        print("warning: 2015 seems to be missing.")

    da = mrms.rainrate.sel(time=v_time)
    da = convert_mm_per_hour_to_in_per_event(da)
    da = da.mean("time")
    
    da.plot.pcolormesh(x="longitude", y="latitude", ax=ax, 
                        cbar_kwargs=dict(pad=0.02, shrink=0.5, label="Event Total Rainfall (inches)"),
                         vmin=0, vmax=np.ceil(float(cbar_high)), cmap='gist_rainbow')
    
    gds_gages_trns_with_rainfall.plot(ax=ax, zorder=2, column="gage_precip_in",
                                     vmin=0, vmax=np.ceil(float(cbar_high)),
                                     cmap='gist_rainbow', edgecolor="black",
                                     missing_kwds=dict(color="none", edgecolor="lightgrey", label = "missing values"))

    for x, y, label in zip(gds_gages_trns_with_rainfall.geometry.x, gds_gages_trns_with_rainfall.geometry.y, gds_gages_trns_with_rainfall.MONITORING):
        xytext=(3, 3)
        if label == "MMPS-140":
            xytext=(-40, 5)
        ax.annotate(label, xy=(x, y), xytext=xytext, 
            textcoords="offset points", color="white", fontsize="x-small",
            fontweight='demibold')

    gds_subs_trns.plot(ax=ax, color="grey", edgecolor="none", alpha = 0.5)
    # gds_gages_trns.plot(ax=ax, zorder=2, color="grey", edgecolor="black")
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    gls = ax.gridlines(draw_labels=True)
    gls.top_labels=False   # suppress top labels
    gls.right_labels=False # suppress right labels
    gds_coast_trns.plot(ax=ax, color='black', zorder=1)
    # ax.spines.right.set_visible(False)
    ax.set_title(title)
    plt.savefig(fname)

def plot_histogram_of_mrms_vs_gage_data(df, title, fname, bins=20, txt_x = 0.02, txt_y = 0.8):
    dif = df.loc[:, 'mrms_precip_in'] - df.loc[:, 'gage_precip_in']
    mean_dif = np.mean(dif)
    interval = st.norm.interval(confidence=0.95, loc=mean_dif, scale=st.sem(dif))
    str_int = str(np.round(interval, 4))
    CI_range = round(abs(interval[0] - interval[1])/2, 5)
    fig, ax = plt.subplots(dpi=300)
    dif.plot.hist(bins=bins, ax=ax)
    plt.text(txt_x, txt_y, ('$\mu=' + str(round(mean_dif, 5)) + '\pm ' + str(CI_range) + r' \ (\alpha=95 \%) $'), transform=ax.transAxes)
    ax.set_title(title)
    ax.set_xlabel("Difference (in)")
    ax.set_ylabel("Count")
    plt.tight_layout()
    plt.savefig(fldr_out_plots+fname)

#%% try plotting all the events
# in order of largest to smallest events
event_max_time = df_gage_and_mrms.loc[:,["event_id", "time"]].groupby(["event_id"]).max()
event_min_time = df_gage_and_mrms.loc[:,["event_id", "time"]].groupby(["event_id"]).min()
event_duration = event_max_time - event_min_time

# removing invalid readings from df_gage_and_mrms
event_timeseries = df_gage_and_mrms.set_index(["event_id", "gage_id", "time"])
m = event_timeseries < 10
event_timeseries = event_timeseries.where(m, np.nan)

# event totals
# a min_count of 4 means there has to be at least an hour's worth of non-NA values
event_totals = event_timeseries.groupby(["event_id", "gage_id"]).sum(min_count=4)

# plot difs at each timestep
event_timeseries_no_na = event_timeseries.dropna()
title = "MRMS minus Gage Depth at Each Timestep \n Where Neither are NA"
fname = "b_mrms_minus_gage_time_series_histogram.png"
plot_histogram_of_mrms_vs_gage_data(event_timeseries_no_na, title, fname)

# plot difs at each timestep with where at least one is non-zero
event_timeseries_no_na_no_zeros = event_timeseries[(event_timeseries.mrms_precip_in>0) | (event_timeseries.gage_precip_in>0)]
event_timeseries_no_na_no_zeros = event_timeseries_no_na_no_zeros.dropna()
title = "MRMS minus Gage Depth at Each Timestep \n Where Neither are NA and At Least One Is Non-zero"
fname = "b_mrms_minus_gage_no_zeros_time_series_histogram.png"
plot_histogram_of_mrms_vs_gage_data(event_timeseries_no_na_no_zeros, title, fname)

# plot differences of event totals
event_totals_no_na = event_totals.dropna()
title = "MRMS minus Gage Event Totals \n Where Neither are NA"
fname = "b_mrms_minus_gage_event_totals_histogram.png"
plot_histogram_of_mrms_vs_gage_data(event_totals_no_na, title, fname, txt_x = 0.5, txt_y = 0.8)


# plot the time series for the top ~30 mrms and gage data
for i in np.arange(10, 30):
    n_largest = i
    top_mrms = event_totals.groupby("event_id").mean().nlargest(n_largest, "mrms_precip_in")
    top_gage = event_totals.groupby("event_id").mean().nlargest(n_largest, "gage_precip_in")
    event_ids = np.unique(np.concatenate((top_gage.index.values, top_mrms.index.values)))
    if len(event_ids) >= 30:
        break

gage_ids = np.unique(event_totals.index.get_level_values(1).values)

 # don't want plots showing up on console
 #%% plotting individual events
 # empty and create folders for each gage
for g_id in gage_ids:
    try:
        shutil.rmtree(fldr_out_plots+g_id)
    except:
        pass
    try:
        os.mkdir(fldr_out_plots+g_id)
    except:
        continue
for g_id, e_id  in itertools.product(gage_ids, event_ids):
    # print("{}, {}".format(g_id, e_id))
    fig, ax = plt.subplots(dpi=150)
    tseries = df_gage_and_mrms[df_gage_and_mrms["event_id"]==e_id].time
    str_start_datetime = ":".join(str(min(tseries)).split(":")[0:2])
    str_end_datetime = ":".join(str(max(tseries)).split(":")[0:2])
    # plt.ioff()
    df_plt = df_gage_and_mrms[df_gage_and_mrms["gage_id"]==g_id]
    df_plt = df_plt[df_plt["event_id"]==e_id]
    df_plt = df_plt.set_index("time")
    df_plt = df_plt.loc[:, ["mrms_precip_in", "gage_precip_in"]]
    # skip when there is no overlapping data
    if (df_plt.sum()[0] == 0) or (df_plt.sum()[1] == 0):
        continue

    df_plt.plot(ax=ax, xlabel="", ylabel="rainfall depth (in)",
                title="gage: {}, event: {} ({} to {})".format(g_id, e_id, str_start_datetime, str_end_datetime))
    plt.tight_layout()
    plt.savefig(fldr_out_plots+g_id + "/e_id_{}".format(e_id))
    plt.close()

#%% plot the mrms data
cbar_high_in = np.ceil(event_totals.mrms_precip_in.max())

try:
    shutil.rmtree(fldr_out_plots+"mrms")
except:
    pass
try:
    os.mkdir(fldr_out_plots+"mrms")
except:
    pass

problem_events = []

for e_id in event_ids:
    tseries = df_gage_and_mrms[df_gage_and_mrms["event_id"]==e_id].time
    str_start_datetime = ":".join(str(min(tseries)).split(":")[0:2])
    str_end_datetime = ":".join(str(max(tseries)).split(":")[0:2])
    gage_totals = event_totals.gage_precip_in.loc[(e_id,)]

    title = "Event {}: {} to {}".format(e_id, str_start_datetime, str_end_datetime)
    fname = fldr_out_plots + "mrms/e_id_{}.png".format(e_id)
    try:
        plot_day_of_mrms_data(gage_totals, tseries, cbar_high_in, title, fname)
    except:
        problem_events.append(title)

print("########################################")
print("Events without corresponding MRMS plot:")
for p in problem_events:
    print(p)
#%% compare with the datasets that were processed the previous way
# df_comp_old = pd.read_csv(fl_csv_prior_processing, index_col = 0, parse_dates=["time"])
