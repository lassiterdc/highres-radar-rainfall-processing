#%% directories 
fldr_stageIV_data = "D:/GDrive/grad_school/StageIV_rainfall/{}/{}.nc" # (year, file name)
fldr_repo = "D:/Dropbox/_GradSchool/_norfolk/highres-radar-rainfall-processing/"
fldr_data = fldr_repo + "data/"
f_ncs_fullres = fldr_data + "mrms_nc_preciprate_fullres_yearlyfiles_atgages/*.nc"
f_ncs_fullres_stageiv = fldr_data + "stage_iv_nc_preciprate_fullres_yearlyfiles_atgages/*.nc"
f_nc_yearly = fldr_data + "mrms_nc_preciprate_yearly_singlefile.nc"

f_nc_atgages = fldr_data + "mrms_nc_preciprate_fullres_atgages.nc"
f_nc_stageiv_at_gages = fldr_data + "stage_iv_nc_preciprate_fullres_yearlyfiles_atgages/"
f_nc_24h_atgages = fldr_data + "mrms_nc_preciprate_24h_atgages.nc"
f_nc_1M_atgages = fldr_data + "mrms_nc_preciprate_monthly_atgages.nc"

f_mrms_csvs_fullres = fldr_data + "mrms_csv_preciprate_fullres_yearlyfiles_atgages/*.csv"
f_stageiv_csvs_fullres = fldr_data + "stage_iv_csv_preciprate_fullres_yearlyfiles_atgages/*.csv"
fl_events = fldr_data + "gage_event_data/rainfall_events_all.csv"
f_csv_mrms_and_gage_events = fldr_data + "/gage_and_mrms_csv_event_preciprate.csv"
f_csv_stageiv_and_gage_events = fldr_data + "/gage_and_stageiv_csv_event_preciprate.csv"


# f_in_nc_24h = "D:/mrms_processing/out_netcdfs/a1_mrms_at_norfolk_24h.nc"
# f_in_nc = "D:/mrms_processing/out_netcdfs/a1_mrms_at_norfolk.nc"
# f_in_csv = "D:/mrms_processing/out_csvs/a_mrms_and_hrsd_event_data.csv"
fldr_out_plots = "D:/Dropbox/_GradSchool/_norfolk/highres-radar-rainfall-processing/plots/"
fldr_out_plots_b = fldr_out_plots + "b_visualizations_of_gage_vs_mrms_data/"
fldr_out_plots_d = fldr_out_plots + "d_visualizations_of_gage_vs_mrms_and_stageIV_data/"
fldr_out_plots_h = fldr_out_plots + "h_annual_statistics/"
f_shp_coast = fldr_data + "geospatial/composite_shoreline_clipped.shp"
f_shp_subcatchments = fldr_data + "geospatial/subcatchments.shp"
f_shp_gages = fldr_data + "geospatial/rain_gages.shp"
f_shp_states = fldr_data + "geospatial/States_shapefile.shp"
f_shp_nexrad_boundary = fldr_repo + "arcpro/shapefiles/nexrad_boundary.shp"

# downloading nexrad station information
fldr_nexrad_data = fldr_data + "nexrad_data/"
f_stations = fldr_nexrad_data + "nexrad_station_metadata.csv"
f_station_data_with_dates = fldr_nexrad_data + "nexrad_station_metadata_with_dates.csv"
url_stations = "https://www.ncei.noaa.gov/access/homr/services/station/{}?date=all" # station id goes here e.g., 30001884
f_nexrad_station_shp = fldr_nexrad_data + "nexrad_stations.shp"



# filepaths for _work scripts
## work_a
fld_mrms_hourly = fldr_data + "mrms_nc_preciprate_hourly_dailyfiles/"
f_stormcat_mrms_hrly_2003 = "D:/Dropbox/_GradSchool/_norfolk/stormy/stochastic_storm_transposition/norfolk/norfolk_mrms_sst_subset_short_2003_hourly_rivanna.nc"
## work b
f_nc_yearly_stageIV = fldr_data + "stageiv_nc_preciprate_yearly_singlefile.nc"


# chunking parameters
size_of_float32 = 32 # bits
MB_per_bit = 1.25e-7
num_lats = 3500
num_lons = 7000
da_chnk_sz = "10000MB" # script seems to succeed when this is 1/4 of the memory allocated (I didn't try to push it)
db_chnk_sz = "5000MB"
dc_chnk_sz = "5000MB" 
ha_chnk_sz = "5000MB" 
hb_chnk_sz = "5000MB" 
i_chnk_sz = "10000MB"

# percentile for colorbar for plotting mrms data
cbar_percentile = 0.98
nearest_int_for_rounding = 2 # inches of rainfall

plt_hb_width = 14
plt_hb_height = plt_hb_width / 1.5

exclude_2013to2014_from_mean_for_anamolies_plot = True

coarsen_to = 0.1 # degree grid

mm_per_in = 25.4
#%% filepaths
def return_colorbar_percentile_for_plotting_gridded_precip_data():
    return cbar_percentile, nearest_int_for_rounding


def return_chunking_parameters(script_prefix, num_lats=num_lats, num_lons=num_lons):
    if script_prefix == "da":
        chnk_sz = da_chnk_sz
    elif script_prefix == "db":
        chnk_sz = db_chnk_sz
    elif script_prefix == "dc":
        chnk_sz = dc_chnk_sz
    elif script_prefix == "ha":
        chnk_sz = ha_chnk_sz
    elif script_prefix == "ha2":
        chnk_sz = ha_chnk_sz
        num_lats = 900
        num_lons = 2100
    elif script_prefix == "hb":
        chnk_sz = hb_chnk_sz
    elif script_prefix == "i":
        chnk_sz = i_chnk_sz
    return chnk_sz, size_of_float32, MB_per_bit, num_lats, num_lons

def return_a1_filepaths():
    return f_ncs_fullres, f_nc_atgages, f_nc_24h_atgages, f_nc_1M_atgages

def return_a2_filepaths():
    return f_mrms_csvs_fullres, f_stageiv_csvs_fullres, fl_events, f_csv_mrms_and_gage_events, f_csv_stageiv_and_gage_events

def return_b_filepaths():
    return f_nc_24h_atgages, f_nc_atgages, f_csv_mrms_and_gage_events, fldr_out_plots_b, f_shp_coast, f_shp_subcatchments, f_shp_gages

def return_c_filepaths():
    return url_stations, fldr_nexrad_data, f_station_data_with_dates, f_stations, f_nexrad_station_shp

def return_d_filepaths():
    return f_nc_24h_atgages, f_nc_atgages, f_nc_stageiv_at_gages, fldr_stageIV_data, f_csv_mrms_and_gage_events, f_csv_stageiv_and_gage_events, fldr_out_plots_d, f_shp_coast, f_shp_subcatchments, f_shp_gages

def return_worka_filepaths():
    return fld_mrms_hourly, f_stormcat_mrms_hrly_2003

def return_workb_filepaths():
    return mm_per_in, coarsen_to, f_nc_yearly, f_nc_yearly_stageIV, f_shp_coast, f_shp_states, f_shp_nexrad_boundary, fldr_out_plots_d

def return_hb_filepaths():
    return mm_per_in, coarsen_to, f_nc_yearly, f_shp_states, f_shp_nexrad_boundary, fldr_out_plots_h
