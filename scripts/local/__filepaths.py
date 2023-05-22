#%% directories 
fldr_stageIV_data = "D:/GDrive/grad_school/StageIV_rainfall/{}/{}.nc" # (year, file name)
fldr_data = "D:/Dropbox/_GradSchool/_norfolk/highres-radar-rainfall-processing/data/"
f_ncs_fullres = fldr_data + "mrms_nc_preciprate_fullres_yearlyfiles_atgages/*.nc"
f_nc_yearly = fldr_data + "mrms_nc_preciprate_yearly_singlefile.nc"

f_nc_atgages = fldr_data + "mrms_nc_preciprate_fullres_atgages.nc"
f_nc_24h_atgages = fldr_data + "mrms_nc_preciprate_24h_atgages.nc"
f_nc_1M_atgages = fldr_data + "mrms_nc_preciprate_monthly_atgages.nc"

f_csvs_fullres = fldr_data + "mrms_csv_preciprate_fullres_yearlyfiles_atgages/*.csv"
fl_events = fldr_data + "gage_event_data/rainfall_events_all.csv"
f_csv_mrms_and_gage_events = fldr_data + "/gage_and_mrms_csv_event_preciprate.csv"


# f_in_nc_24h = "D:/mrms_processing/out_netcdfs/a1_mrms_at_norfolk_24h.nc"
# f_in_nc = "D:/mrms_processing/out_netcdfs/a1_mrms_at_norfolk.nc"
# f_in_csv = "D:/mrms_processing/out_csvs/a_mrms_and_hrsd_event_data.csv"
fldr_out_plots = "D:/Dropbox/_GradSchool/_norfolk/highres-radar-rainfall-processing/plots/b_visualizations_of_gage_vs_mrms_data/"
f_shp_coast = fldr_data + "geospatial/composite_shoreline_clipped.shp"
f_shp_subcatchments = fldr_data + "geospatial/subcatchments.shp"
f_shp_gages = fldr_data + "geospatial/rain_gages.shp"


# downloading nexrad station information
fldr_nexrad_data = fldr_data + "nexrad_data/"
f_stations = fldr_nexrad_data + "nexrad_station_metadata.csv"
f_station_data_with_dates = fldr_nexrad_data + "nexrad_station_metadata_with_dates.csv"
url_stations = "https://www.ncei.noaa.gov/access/homr/services/station/{}?date=all" # station id goes here e.g., 30001884
f_nexrad_station_shp = fldr_nexrad_data + "nexrad_stations.shp"

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

def return_c_filepaths():
    return url_stations, fldr_nexrad_data, f_station_data_with_dates, f_stations, f_nexrad_station_shp

def return_worka_filepaths():
    return fld_mrms_hourly, f_stormcat_mrms_hrly_2003

def return_workb_filepaths():
    return fldr_stageIV_data, f_nc_yearly
