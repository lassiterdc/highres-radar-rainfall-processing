
# inputs
fldr_data = "D:/Dropbox/_GradSchool/_norfolk/highres-radar-rainfall-processing/data/"
f_ncs_fullres = fldr_data + "mrms_nc_preciprate_fullres_yearlyfiles_atgages/*.nc"

# outputs
f_nc_atgages = fldr_data + "mrms_nc_preciprate_fullres_atgages.nc"
f_nc_24h_atgages = fldr_data + "mrms_nc_preciprate_24h_atgages.nc"
f_nc_1M_atgages = fldr_data + "mrms_nc_preciprate_monthly_atgages.nc"

# a2 specific files
f_csvs_fullres = fldr_data + "mrms_csv_preciprate_fullres_yearlyfiles_atgages/*.csv"
fl_events = fldr_data + "gage_event_data/rainfall_events_all.csv"
f_csv_mrms_and_gage_events = fldr_data + "/gage_and_mrms_csv_event_preciprate.csv"

# b specific files
# f_in_nc_24h = "D:/mrms_processing/out_netcdfs/a1_mrms_at_norfolk_24h.nc"
# f_in_nc = "D:/mrms_processing/out_netcdfs/a1_mrms_at_norfolk.nc"
# f_in_csv = "D:/mrms_processing/out_csvs/a_mrms_and_hrsd_event_data.csv"
fldr_out_plots = "D:/Dropbox/_GradSchool/_norfolk/highres-radar-rainfall-processing/plots/b_visualizations_of_gage_vs_mrms_data/"
f_shp_coast = fldr_data + "geospatial/composite_shoreline_clipped.shp"
f_shp_subcatchments = fldr_data + "geospatial/subcatchments.shp"
f_shp_gages = fldr_data + "geospatial/rain_gages.shp"

# filepaths for _work scripts
## work_a
fld_mrms_hourly = "D:/Dropbox/_GradSchool/_norfolk/highres-radar-rainfall-processing/data/mrms_nc_preciprate_hourly_dailyfiles/"
f_stormcat_mrms_hrly_2003 = "D:/Dropbox/_GradSchool/_norfolk/stormy/stochastic_storm_transposition/norfolk/norfolk_mrms_sst_subset_short_2003_hourly_rivanna.nc"



#%% filepaths for script a1
def return_a1_filepaths():
    return f_ncs_fullres, f_nc_atgages, f_nc_24h_atgages, f_nc_1M_atgages

def return_a2_filepaths():
    return f_csvs_fullres, fl_events, f_csv_mrms_and_gage_events

def return_b_filepaths():
    return f_nc_24h_atgages, f_nc_atgages, f_csv_mrms_and_gage_events, fldr_out_plots, f_shp_coast, f_shp_subcatchments, f_shp_gages

def return_worka_filepaths():
    return fld_mrms_hourly, f_stormcat_mrms_hrly_2003

