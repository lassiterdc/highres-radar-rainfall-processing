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
from glob import glob
from __filepaths import return_d_filepaths
import sys
import matplotlib
import matplotlib.patheffects as PathEffects
import math
# Set matplotlib parameters
plt.rcParams['text.usetex'] = True
# default text size
plt.rc('font', size=11)
# Set the axes title font size
plt.rc('axes', titlesize=16)
# Set the font size of the figure title
plt.rc('figure', titlesize=16)
# Set the axes labels font size
plt.rc('axes', labelsize=14)
# Set the font size for x tick labels
plt.rc('xtick', labelsize=12)
# Set the font size for y tick labels
plt.rc('ytick', labelsize=12)
# Set the legend font size
plt.rc('legend', fontsize=16)



figwidth, figheight = (7.5, 6)

n_largest = 15 # largest events to inspect

exclude_2012_2013_and_2014 = True

inches_per_mm = 1/25.4

lst_event_comparisons = [107, 5]
lst_gage_comparisons = ['MMPS-149', 'MMPS-163']
# f_in_nc_24h = "D:/mrms_processing/out_netcdfs/a1_mrms_at_norfolk_24h.nc"
# f_in_nc = "D:/mrms_processing/out_netcdfs/a1_mrms_at_norfolk.nc"
# f_in_csv = "D:/mrms_processing/out_csvs/a_mrms_and_hrsd_event_data.csv"
# fldr_out_plots = "D:/mrms_processing/plots/b_visualizations_of_gage_vs_mrms_data/"
# f_shp_coast = "D:/mrms_processing/data/shapefiles/composite_shoreline_clipped.shp"
# f_shp_subcatchments = "D:/mrms_processing/data/shapefiles/subcatchments.shp"
# f_shp_gages = 'D:/mrms_processing/data/gage_hrsd/rain_gages.shp'

f_nc_24h_atgages, f_nc_atgages, fldr_nc_stageiv_at_gages, fldr_stageIV_data, f_csv_mrms_and_gage_events, f_csv_stageiv_and_gage_events, fldr_out_plots, f_shp_coast, f_shp_subcatchments, f_shp_gages = return_d_filepaths()
#%% analyze event totals
# loading and preprocessing data
df_gage_and_mrms = pd.read_csv(f_csv_mrms_and_gage_events, parse_dates=["time"], infer_datetime_format=True, index_col=0)
df_gage_and_stageiv = pd.read_csv(f_csv_stageiv_and_gage_events, parse_dates=["time"], infer_datetime_format=True, index_col=0)

# rename precip columns
df_gage_and_mrms.rename(columns = dict(precip_in = "mrms_precip_in"), inplace=True)
df_gage_and_stageiv.rename(columns = dict(precip_in = "stageiv_precip_in"), inplace=True)
# join gage and mrms data into single dataframe
df_mrms_stageiv_gage = pd.merge(df_gage_and_stageiv, df_gage_and_mrms, on = ["time", "gage_id", "event_id", "gage_precip_in"], how = "inner")

# processing event data
# in order of largest to smallest events
event_max_time = df_mrms_stageiv_gage.loc[:,["event_id", "time"]].groupby(["event_id"]).max()
event_min_time = df_mrms_stageiv_gage.loc[:,["event_id", "time"]].groupby(["event_id"]).min()
event_duration = event_max_time - event_min_time

# removing invalid readings from df_mrms_stageiv_gage
event_timeseries = df_mrms_stageiv_gage.set_index(["event_id", "gage_id", "time"])
m = event_timeseries < 10
event_timeseries = event_timeseries.where(m, np.nan)

# event totals
# a min_count of 4 means there has to be at least an hour's worth of non-NA values
event_totals = event_timeseries.groupby(["event_id", "gage_id"]).sum(min_count=4)

event_timeseries_no_na_no_zeros = event_timeseries[(event_timeseries.mrms_precip_in>0) | (event_timeseries.gage_precip_in>0)]
event_timeseries_no_na_no_zeros = event_timeseries_no_na_no_zeros.dropna()

event_totals_no_na = event_totals.dropna()

top_mrms = event_totals.groupby("event_id").mean().nlargest(n_largest, "mrms_precip_in")
event_ids = top_mrms.index.values
# event_mrms_ranks = np.arange(1, n_largest+1)
gage_ids = np.unique(event_totals.index.get_level_values(1).values)

# loading netcdf data
def create_time_mask_to_exclude_target_years(tseries):
    dti_years = pd.DatetimeIndex(tseries).year
    exclude = dti_years.where((dti_years != 2012) & (dti_years != 2013) & (dti_years != 2014))
    mask = exclude >= 0 # this yields a False at nan locations
    return mask

mrms_1d = xr.open_dataset(f_nc_24h_atgages).load()
mrms = xr.open_dataset(f_nc_atgages, chunks={"time":"1000MB"})

if exclude_2012_2013_and_2014 == True:
    mask = create_time_mask_to_exclude_target_years(df_gage_and_mrms.time)
    df_gage_and_mrms = df_gage_and_mrms[mask]

    mask = create_time_mask_to_exclude_target_years(mrms_1d.time.values)
    times_idx = (np.where(mask == True))[0]
    mrms_1d = mrms_1d.isel(time=times_idx)

    mask = create_time_mask_to_exclude_target_years(mrms.time.values)
    times_idx = (np.where(mask == True))[0]
    mrms = mrms.isel(time=times_idx)

lst_yrs = []
lst_f_stageiv = []
for yr in df_gage_and_mrms.time.dt.year.unique():
    if exclude_2012_2013_and_2014 == True:
        if yr in [2012, 2013, 2014]:
            continue
        else:
            lst_yrs.append(yr)
    for f in glob(fldr_nc_stageiv_at_gages + "{}.nc".format(yr)):
        lst_f_stageiv.append(f)

stage_iv =  xr.open_mfdataset(lst_f_stageiv, engine='netcdf4')

# stage_iv['outlon'] =  stage_iv.longitude.isel(time=0)+360
# stage_iv['outlat'] =  stage_iv.latitude.isel(time=0)
# stage_iv = stage_iv.drop_vars(["latitude", "longitude"])
# stage_iv = stage_iv.rename_dims(dims_dict=dict(outlat = "latitude", outlon="longitude"))
# stage_iv = stage_iv.rename(dict(outlat = "latitude", outlon = "longitude"))

# load shapefiles
gds_subs = gpd.read_file(f_shp_subcatchments)
gds_gages = gpd.read_file(f_shp_gages)
gds_coast = gpd.read_file(f_shp_coast)

proj = ccrs.PlateCarree()
gds_subs_trns = gds_subs.to_crs(proj)

gds_gages_trns = gds_gages.to_crs(proj)
gds_coast_trns = gds_coast.to_crs(proj)

# make sure plotting folder exists
# try:
#     shutil.rmtree(fldr_out_plots)
# except:
#     pass
try:
    os.mkdir(fldr_out_plots)
except:
    pass

#%% define plotting functions
def convert_mm_per_hour_to_in_per_event(da):
    # tseries = da.time.values
    # https://pandas.pydata.org/pandas-docs/stable/user_guide/timeseries.html#offset-aliases
    tseries = pd.DatetimeIndex(da.time.values)
    freq = tseries.inferred_freq
    period = 1
    try:
        period=int(freq.split("T")[0])
        freq="T"
    except:
        pass
    if freq is not None:
        hrs_in_event = pd.Timedelta(len(tseries), freq)*period / np.timedelta64(1, 'h')
    # hrs_in_event = (max(tseries) - min(tseries)) / np.timedelta64(1, 'h')
    da = da * hrs_in_event * inches_per_mm # mm/hr * hrs in event * inches_per_mm
    return da

def plot_day_of_mrms_data(gage_totals, tseries, cbar_high, title, fname, gds_subs_trns=gds_subs_trns, gds_gages_trns=gds_gages_trns, gds_coast_trns=gds_coast_trns):
    fig, ax = plt.subplots(figsize = [figwidth, figheight], subplot_kw=dict(projection=proj), dpi=500)
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
                         vmin=0, vmax=np.ceil(float(cbar_high)), cmap=cmap)
    
    gds_gages_trns_with_rainfall.plot(ax=ax, zorder=2, column="gage_precip_in",
                                     vmin=0, vmax=np.ceil(float(cbar_high)),
                                     cmap=cmap, edgecolor="black",
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
    return ax, xlim, ylim

def plot_compare_stageIV_and_mrms(gage_totals, tseries, cbar_high, title, fname, cmap = "plasma", xlim=(-76.38, -76.12), ylim=(36.79, 36.97), gds_subs_trns=gds_subs_trns, gds_gages_trns=gds_gages_trns, gds_coast_trns=gds_coast_trns):
    fig, axes = plt.subplots(1,2, sharey="row", sharex = True, figsize = [figwidth*2, figheight*0.8], subplot_kw=dict(projection=proj), dpi=500)
    ######### PLOT MRMS #########
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
    
    da.plot.pcolormesh(x="longitude", y="latitude", ax=axes[0],
                         vmin=0, vmax=np.ceil(float(cbar_high)), cmap=cmap,
                         add_colorbar = False)
    
    gds_gages_trns_with_rainfall.plot(ax=axes[0], zorder=2, column="gage_precip_in",
                                     vmin=0, vmax=np.ceil(float(cbar_high)),
                                     cmap=cmap, edgecolor="black",
                                     missing_kwds=dict(color="none", edgecolor="lightgrey", label = "missing values"))
    
    # setting x and y extents:
    # ax_tmp = gds_gages_trns_with_rainfall.plot(zorder=2, column="gage_precip_in",
    #                                 vmin=0, vmax=np.ceil(float(cbar_high)),
    #                                 cmap=cmap, edgecolor="black",
    #                                 missing_kwds=dict(color="none", edgecolor="lightgrey", label = "missing values"))

    for x, y, label in zip(gds_gages_trns_with_rainfall.geometry.x, gds_gages_trns_with_rainfall.geometry.y, gds_gages_trns_with_rainfall.MONITORING):
        xytext=(3, 6)
        if label == "MMPS-140":
            xytext=(-40, 6)
        txt = axes[0].annotate(label, xy=(x, y), xytext=xytext, 
            textcoords="offset points", color="black", fontsize=11,
            fontweight='heavy')
        txt.set_path_effects([PathEffects.withStroke(linewidth=3.8, foreground='w')])
        # plt.draw()

    gds_subs_trns.plot(ax=axes[0], color="grey", edgecolor="none", alpha = 0.5)
    # gds_gages_trns.plot(ax=ax, zorder=2, color="grey", edgecolor="black")
    # xlim = ax_tmp.get_xlim()
    # ylim = ax_tmp.get_ylim()
    # print(xlim)
    # print(ylim)
    axes[0].set_xlim(xlim)
    axes[0].set_ylim(ylim)
    gls = axes[0].gridlines(draw_labels=True)
    gls.top_labels=False   # suppress top labels
    gls.right_labels=False # suppress right labels
    gds_coast_trns.plot(ax=axes[0], color='black', zorder=1)
    # ax.spines.right.set_visible(False)
    axes[0].set_title("MRMS")

    ######### PLOT STAGE IV #########
    # modifying gages geodatabase
    gage_totals_renamed = gage_totals.reset_index().rename(columns={"gage_id":"MONITORING"})
    gds_gages_trns_with_rainfall = gds_gages_trns.merge(gage_totals_renamed, on="MONITORING")
    
    # subsetting based on time
    tseries_min = min(tseries.values)
    tseries_max = max(tseries.values)
    # v_time = mrms.rainrate.time.values
    v_time = stage_iv.rainrate.time.values
    v_time = v_time[(v_time>=tseries_min) & (v_time<=tseries_max)]

    # troubleshooting 2015 not showing up 
    # dti = pd.DatetimeIndex(mrms.rainrate.time.values)
    dti = pd.DatetimeIndex(stage_iv.rainrate.time.values)
    if sum(dti.year == 2015) == 0:
        print("warning: 2015 seems to be missing.")

    # da = mrms.rainrate.sel(time=v_time)
    da = stage_iv.rainrate.sel(time=v_time)
    # WORK
    # da = convert_mm_per_hour_to_in_per_event(da)
    # da = da.mean("time")
    # END WORK
    tseries = pd.DatetimeIndex(da.time.values)
    freq = tseries.inferred_freq
    if freq is not None:
        hrs_in_event = pd.Timedelta(len(tseries), freq) / np.timedelta64(1, 'h')
    else: 
        import sys
        sys.exit("frequency could not be autodetected")
    da_mean_intensity = da.mean("time") # mm/hr
    da_event_total_inches = da_mean_intensity *  hrs_in_event # mm/hr * event duration in hours = mm
    da_event_total_inches = da_event_total_inches * inches_per_mm # mm * inches_per_mm = in
    
    # WORK
    # da_event_total_inches.plot.pcolormesh(x="longitude", y="latitude", ax=axes[1], 
    #                     cbar_kwargs=dict(pad=0.02, shrink=0.5, label="Event Total Rainfall (inches)"),
    #                      vmin=0, vmax=np.ceil(float(cbar_high)), cmap=cmap)
    
    im = da_event_total_inches.plot.pcolormesh(x="longitude", y="latitude", ax=axes[1], 
                        vmin=0, vmax=np.ceil(float(cbar_high)), cmap=cmap, add_colorbar = False)
    
    cbar_ax = fig.add_axes([0.92, 0.13, 0.01, 0.7])
    fig.colorbar(im, cax=cbar_ax, pad=0.02, shrink=0.5, label="Event Total Rainfall (inches)")
    # END WORK
    gds_gages_trns_with_rainfall.plot(ax=axes[1], zorder=2, column="gage_precip_in",
                                     vmin=0, vmax=np.ceil(float(cbar_high)),
                                     cmap=cmap, edgecolor="black",
                                     missing_kwds=dict(color="none", edgecolor="lightgrey", label = "missing values"))

    for x, y, label in zip(gds_gages_trns_with_rainfall.geometry.x, gds_gages_trns_with_rainfall.geometry.y, gds_gages_trns_with_rainfall.MONITORING):
        xytext=(3, 6)
        if label == "MMPS-140":
            xytext=(-40, 6)
        txt = axes[1].annotate(label, xy=(x, y), xytext=xytext, 
            textcoords="offset points", color="black", fontsize=11,
            fontweight='heavy')
        txt.set_path_effects([PathEffects.withStroke(linewidth=3.8, foreground='w')])



    gds_subs_trns.plot(ax=axes[1], color="grey", edgecolor="none", alpha = 0.5)
    # gds_gages_trns.plot(ax=ax, zorder=2, color="grey", edgecolor="black")
    # xlim = ax.get_xlim()
    # ylim = ax.get_ylim()
    
    axes[1].set_xlim(xlim)
    axes[1].set_ylim(ylim)
    gls = axes[1].gridlines(draw_labels=True)
    gls.top_labels=False   # suppress top labels
    gls.right_labels=False # suppress right labels
    gds_coast_trns.plot(ax=axes[1], color='black', zorder=1)
    # ax.spines.right.set_visible(False)
    axes[1].set_title("Stage IV")
    fig.suptitle(title)
    # labs = axes[1].axes.get_yticklabels()
    # vals = axes[1].axes.get_yticks()
    # print(labs)
    # new_labs = []
    # for lab in labs:
    #     new_labs.append("")
    # axes[1].axes.set_yticks(vals)
    # axes[1].axes.set_yticklabels(new_labs)
    # print(vals)
    # print(new_labs)
    # axes[1].set_ylabel("")
    # axes[1].axes.set_yticks([])
    # plt.tight_layout()
    plt.savefig(fname)

# def plot_histogram_of_mrms_vs_gage_data(df, title, fname, bins=20, txt_x = 0.02, txt_y = 0.8):
#     dif = df.loc[:, 'mrms_precip_in'] - df.loc[:, 'gage_precip_in']
#     mean_dif = np.mean(dif)
#     interval = st.norm.interval(confidence=0.95, loc=mean_dif, scale=st.sem(dif))
#     str_int = str(np.round(interval, 4))
#     CI_range = round(abs(interval[0] - interval[1])/2, 5)
#     fig, ax = plt.subplots(dpi=300)
#     dif.plot.hist(bins=bins, ax=ax)
#     plt.text(txt_x, txt_y, ('$\mu=' + str(round(mean_dif, 5)) + '\pm ' + str(CI_range) + r' \ (\alpha=95 \%) $'), transform=ax.transAxes)
#     ax.set_title(title)
#     ax.set_xlabel("Difference (in)")
#     ax.set_ylabel("Count")
#     plt.tight_layout()
#     plt.savefig(fldr_out_plots+fname)

# def compute_ci(mean_dif, dif):
#     interval = st.norm.interval(confidence=0.95, loc=mean_dif, scale=st.sem(dif))
#     str_int = str(np.round(interval, 4))
#     CI_range = round(abs(interval[0] - interval[1])/2, 5)
#     return CI_range

def compute_ci_bootstrap(dif, alpha, round_to):
    from scipy.stats import bootstrap
    data = (dif.values,)  # samples must be in a sequence
    res = bootstrap(data, np.mean, confidence_level=alpha)
    CI_low = round(res.confidence_interval.low, round_to)
    CI_high = round(res.confidence_interval.high, round_to)
    return CI_low, CI_high

def plot_histogram_of_mrms_vs_stageiv_vs_gage_data(df, title, fname, quant = None, bin_width = 0.2, txt_x = 0.02, txt_y = 0.8):
    dif_mrms = df.loc[:, 'mrms_precip_in'] - df.loc[:, 'gage_precip_in']
    dif_stageiv = df.loc[:, 'stageiv_precip_in'] - df.loc[:, 'gage_precip_in']
    mean_dif_mrms = np.mean(dif_mrms)
    mean_dif_stageiv = np.mean(dif_stageiv)

    CI_low_mrms, CI_high_mrms = compute_ci_bootstrap(dif_mrms, alpha = 0.95, round_to = 4)
    CI_low_stageiv, CI_high_stageiv = compute_ci_bootstrap(dif_stageiv, alpha = 0.95, round_to = 4)

    # define bins so that one is centered on zero and so it contains all the stageiv and mrms data
    all_difs = np.concatenate((dif_mrms, dif_stageiv))

    if quant is not None:
        q_low = (1 - quant) / 2
        q_high = 1 -    q_low
        min_bin = np.quantile(all_difs, q_low)
        max_bin = np.quantile(all_difs, q_high)
    else:
        min_bin = np.floor(min(all_difs) * 100)/100.0 # rounding to the nearest 2 decimals
        max_bin = np.ceil(max(all_difs) * 100)/100.0
    bin_left_edge_next_to_0 = 0 - bin_width/2
    bin_left_edge_right_of_0 = 0 + bin_width/2
    left_bins = np.arange(start = abs(bin_left_edge_next_to_0), stop = abs(min_bin), step = bin_width) * -1
    right_bins = np.arange(start = bin_left_edge_right_of_0, stop = max_bin+bin_width, step = bin_width)

    bin_edges = np.concatenate((np.sort(left_bins), right_bins))


    fig, ax = plt.subplots(dpi=500)
    hist_outs_mrms = ax.hist(x=dif_mrms, bins=bin_edges, color = 'blue', alpha = 0.5, label = 'MRMS')
    hist_outs_stageiv = ax.hist(x=dif_stageiv, bins=bin_edges, color = 'green', alpha = 0.5, label = "Stage IV")

    # plt.text(txt_x, txt_y, ('$\mu_{MRMS}=' + str(round(mean_dif_mrms, 5)) + '\pm ' + str(CI_range_mrms) + r' \ (\alpha=95 \%) $'), transform=ax.transAxes)
    # plt.text(txt_x, txt_y-0.07, ('$\mu_{StageIV}=' + str(round(mean_dif_stageiv, 5)) + '\pm ' + str(CI_range_stageiv) + r' \ (\alpha=95 \%) $'), transform=ax.transAxes)

    plt.text(txt_x, txt_y, ('$' + str(CI_low_mrms) + ' \le ' +'\mu_{MRMS}'+ ' \le ' + str(CI_high_mrms) + r' \ (\alpha=95 \%) $'), transform=ax.transAxes)
    plt.text(txt_x, txt_y-0.07, ('$' + str(CI_low_stageiv) + ' \le ' +'\mu_{Stage IV}'+ ' \le ' + str(CI_high_stageiv) + r' \ (\alpha=95 \%) $'), transform=ax.transAxes)

    ax.set_title(title, horizontalalignment='center',verticalalignment='center')
    ax.set_xlabel("Difference (in)")
    ax.set_ylabel("Count")
    ax.legend()
    plt.tight_layout()
    plt.savefig(fldr_out_plots+fname)

#%% plotting mrms, stageIV, and gage comparison


title = "Gridded Rainfall Minus Gage Rainfall at Non-Zero Timesteps"
quant = 0.98 # this cuts off the tails of the histogram so you can see the patterns in the data easier
fname = "a_mrms_stageiv_minus_gage_no_zeros_time_series_histogram{}.svg".format(str(quant))
plot_histogram_of_mrms_vs_stageiv_vs_gage_data(event_timeseries_no_na_no_zeros, title, fname, bin_width = 0.03, quant = quant, txt_x = 0.02, txt_y = 0.75)

# title = "Difference between MRMS and Stage IV with Gage at Each Timestep \n Where None are NA, at Least One Is Non-zero, and the intensity is at least 1\" per hour"
# fname = "a_mrms_stageiv_minus_gage_no_zeros_time_series_histogram.png"
# plot_histogram_of_mrms_vs_stageiv_vs_gage_data(event_timeseries_no_na_no_zeros, title, fname)


title = "Event Total Gridded Rainfall Minus Gage Rainfall"
quant = 0.99
fname = "a_mrms_stageiv_minus_gage_event_totals_histogram_quant{}.svg".format(str(quant))
plot_histogram_of_mrms_vs_stageiv_vs_gage_data(event_totals_no_na, title, fname, quant = quant, txt_x = 0.49, txt_y = 0.7)
#%% create plots

# plot difs at each timestep with where at least one is non-zero

# title = "MRMS minus Gage Depth at Each Timestep \n Where Neither are NA and At Least One Is Non-zero"
# fname = "b_mrms_minus_gage_no_zeros_time_series_histogram.png"
# plot_histogram_of_mrms_vs_gage_data(event_timeseries_no_na_no_zeros, title, fname)

# # plot differences of event totals

# title = "MRMS minus Gage Event Totals \n Where Neither are NA"
# fname = "b_mrms_minus_gage_event_totals_histogram.png"
# plot_histogram_of_mrms_vs_gage_data(event_totals_no_na, title, fname, txt_x = 0.5, txt_y = 0.8)


# plot the time series for the top 30 mrms data events
# for i in np.arange(10, 30):
#     n_largest = i
#     top_mrms = event_totals.groupby("event_id").mean().nlargest(n_largest, "mrms_precip_in")
#     top_gage = event_totals.groupby("event_id").mean().nlargest(n_largest, "gage_precip_in")
#     event_ids = np.unique(np.concatenate((top_gage.index.values, top_mrms.index.values)))
#     if len(event_ids) >= 30:
#         break


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
# determine plot widths so that 1 hour corresponds to the same horizontal width on a plot
lst_durations = []
for g_id, e_id  in itertools.product(lst_gage_comparisons, lst_event_comparisons):
    event_index = int(np.where(event_ids==e_id)[0])+1
    tseries = df_mrms_stageiv_gage[df_mrms_stageiv_gage["event_id"]==e_id].time
    duration_h = (max(tseries) - min(tseries)) / np.timedelta64(1, 'h')
    lst_durations.append(duration_h)

dur_max = max(lst_durations)
max_fig_width = 15
min_fig_width = 5
fig_width_per_hr = max_fig_width / dur_max
fig_height = 5

# for g_id, e_id  in itertools.product(gage_ids, event_ids):
for g_id, e_id  in itertools.product(lst_gage_comparisons, lst_event_comparisons):
    event_index = int(np.where(event_ids==e_id)[0])+1
    # print("{}, {}".format(g_id, e_id))
    tseries = df_mrms_stageiv_gage[df_mrms_stageiv_gage["event_id"]==e_id].time
    duration_h = (max(tseries) - min(tseries)) / np.timedelta64(1, 'h')
    fig, ax = plt.subplots(dpi=500)
    str_start_datetime = ":".join(str(min(tseries)).split(":")[0:2])
    str_end_datetime = ":".join(str(max(tseries)).split(":")[0:2])
    # plt.ioff()
    df_plt = df_mrms_stageiv_gage[df_mrms_stageiv_gage["gage_id"]==g_id]
    df_plt = df_plt[df_plt["event_id"]==e_id]
    df_plt = df_plt.set_index("time")
    df_plt = df_plt.rename(columns = dict(mrms_precip_in = "MRMS", stageiv_precip_in = "Stage IV", gage_precip_in = "Gage"))
    df_plt = df_plt.loc[:, ["MRMS", "Stage IV", "Gage"]]
    # skip when there is no overlapping data
    if (df_plt.sum()[0] == 0) or (df_plt.sum()[1] == 0):
        continue
    
    fig_width = max(min_fig_width, fig_width_per_hr*duration_h)

    matplotlib.style.use('seaborn-colorblind')
    df_plt.plot(ax=ax, xlabel="", ylabel="Rainfall Depth (in)",
                # title="Gage: {}, Event: {} ({} to {})".format(g_id, e_id, str_start_datetime, str_end_datetime),
                title="Gage: {}, Event: {}".format(g_id, e_id),
                style = ['-', '--', ":"], figsize=(fig_width, fig_height))

    linewidths = 1.4
    for line in ax.get_lines():
        # https://matplotlib.org/stable/api/_as_gen/matplotlib.lines.Line2D.html
        if line.get_label() == 'Gage':
            line.set_linewidth(linewidths*1.2)
            line.set(zorder = 80)
            # line.set_dashes((0, (1, 1)))
            # line.set_color('black')
            line.set_linestyle((0, (1, 1.3)))
        if line.get_label() == 'MRMS':
            line.set_linewidth(linewidths*1.2)
            line.set_color('blue')
            line.set_alpha(0.5)
        if line.get_label() == 'Stage IV':
            line.set_linewidth(linewidths*1.3)
            line.set(zorder = 100)
            line.set_color('green')
            line.set_alpha(0.5)
        ax.legend()

    plt.tight_layout()
    plt.savefig(fldr_out_plots+g_id + "/{}_e_id_{}.svg".format(event_index, e_id))
    # plt.close()

#%% plot the mrms data
# cbar_high_in = np.ceil(event_totals.mrms_precip_in.max())

# try:
#     shutil.rmtree(fldr_out_plots+"mrms")
# except:
#     pass
# try:
#     os.mkdir(fldr_out_plots+"mrms")
# except:
#     pass

# problem_events = []

# for e_id in event_ids:
# # for e_id in event_ids[0:1]:
#     event_index = int(np.where(event_ids==e_id)[0])+1
#     tseries = df_gage_and_mrms[df_gage_and_mrms["event_id"]==e_id].time
#     str_start_datetime = ":".join(str(min(tseries)).split(":")[0:2])
#     str_end_datetime = ":".join(str(max(tseries)).split(":")[0:2])
#     gage_totals = event_totals.gage_precip_in.loc[(e_id,)]

#     title = "Event {}: {} to {}".format(e_id, str_start_datetime, str_end_datetime)
#     fname = fldr_out_plots + "mrms/{}_e_id_{}.png".format(event_index, e_id)
#     # try:
#     #     xlim, ylim = plot_day_of_mrms_data(gage_totals, tseries, cbar_high_in, title, fname)
#     # except:
#     #     problem_events.append(title)

#     xlim, ylim = plot_day_of_mrms_data(gage_totals, tseries, cbar_high_in, title, fname)

# print("########################################")
# print("Events without corresponding MRMS plot:")
# for p in problem_events:
#     print(p)


#%% plot comparison of stage IV and MRMS data
cbar_high_in = np.ceil(event_totals.mrms_precip_in.max())
fldr = fldr_out_plots+"mrms_vs_stageiv/"

try:
    shutil.rmtree(fldr)
except:
    pass
try:
    os.mkdir(fldr)
except:
    pass

problem_events = []

# plt.rc('font', size=12)

lst_event_comparisons = [107]

# for e_id in event_ids:
for e_id in lst_event_comparisons:
# for e_id in event_ids[0:1]:
    event_index = int(np.where(event_ids==e_id)[0])+1
    tseries = df_gage_and_mrms[df_gage_and_mrms["event_id"]==e_id].time
    str_start_datetime = ":".join(str(min(tseries)).split(":")[0:2])
    str_end_datetime = ":".join(str(max(tseries)).split(":")[0:2])
    gage_totals = event_totals.gage_precip_in.loc[(e_id,)]

    title = "Event {}: {} to {}".format(e_id, str_start_datetime, str_end_datetime)
    
    # try:
    #     xlim, ylim = plot_day_of_mrms_data(gage_totals, tseries, cbar_high_in, title, fname)
    # except:
    #     problem_events.append(title)
    # fname = fldr + "plasma_{}_e_id_{}.png".format(event_index, e_id)
    # plot_compare_stageIV_and_mrms(gage_totals, tseries, cbar_high_in, title, fname, cmap = "plasma")

    fname = fldr + "gist_rainbow_{}_e_id_{}.svg".format(event_index, e_id)
    plot_compare_stageIV_and_mrms(gage_totals, tseries, cbar_high_in, title, fname, cmap = "gist_rainbow")
