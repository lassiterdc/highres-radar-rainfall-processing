#!/bin/bash

#%% code examples
# example to call array
# echo ${assar_dirs[working]}
# concatenate
#  echo "${assar_dirs[working]}${assar_dirs[raw_nssl]}"

# example concatenating assar_dirs[hpc_c_py]=${assar_dirs[hpc_scripts]}"_c_running_rainyday_mrms_hourly.py"

# declare associative array
declare -A assar_dirs

# populate associative array with folder and filepaths; this should be the only thing that has to be changed
assar_dirs[repo]="/project/quinnlab/dcl3nd/norfolk/highres-radar-rainfall-processing/"
assar_dirs[stageiv_rainfall]="/project/quinnlab/dcl3nd/norfolk/stormy/data/climate/StageIV_rainfall/"

## paths relative to repo
### hpc script output folder
assar_dirs[hpc_scripts]=${assar_dirs[repo]}"scripts/hpc/"
assar_dirs[hpc_outputs]=${assar_dirs[hpc_scripts]}"_script_outputs/"
assar_dirs[hpc_errors]=${assar_dirs[hpc_scripts]}"_script_errors/"


# data
assar_dirs[hpc_data]=${assar_dirs[repo]}"data/"

### raw data
assar_dirs[raw_nssl]=${assar_dirs[hpc_data]}"raw_data/mrms_grib_nssl/"
assar_dirs[raw_mrms_quantized]=${assar_dirs[hpc_data]}"raw_data/mrms_nc_quant/"
assar_dirs[raw_mrms]=${assar_dirs[hpc_data]}"raw_data/mrms_grib_mesonet/"

### scripts
assar_dirs[hpc_da]=${assar_dirs[hpc_scripts]}"_da_cmbn_to_dly_ncs_frmtd_for_RainyDay.py"
assar_dirs[hpc_db]=${assar_dirs[hpc_scripts]}"_db_resampling_to_hourly_and_daily_timesteps.py"
assar_dirs[hpc_dc]=${assar_dirs[hpc_scripts]}"_dc_combining_daily_totals_in_annual_netcdfs.py"
assar_dirs[hpc_ha]=${assar_dirs[hpc_scripts]}"_ha_generate_annual_statistics_netcdfs.py"
assar_dirs[hpc_ha2]=${assar_dirs[hpc_scripts]}"_ha2_generate_annual_statistics_netcdfs_stageIV.py"
assar_dirs[hpc_hb]=${assar_dirs[hpc_scripts]}"_hb_generate_annual_statistics_plots.py"
assar_dirs[hpc_i]=${assar_dirs[hpc_scripts]}"_i_extract_mrms_at_gages.py"

### scratch folders
assar_dirs[scratch_zarrs]=${assar_dirs[hpc_data]}"_scratch/zarrs/"
assar_dirs[scratch_gribs]=${assar_dirs[hpc_data]}"_scratch/gribs/"

### Outputs:
#### outputs: processed data (format is [source]_[type]_[data]_[timestep]_[file aggregation]_[other qualifiers])
assar_dirs[out_fullres_dailyfiles]=${assar_dirs[hpc_data]}"mrms_nc_preciprate_fullres_dailyfiles/"
assar_dirs[out_hourly_dailyfiles]=${assar_dirs[hpc_data]}"mrms_nc_preciprate_hourly_dailyfiles/"
assar_dirs[out_daily_dailyfiles]=${assar_dirs[hpc_data]}"mrms_nc_preciprate_daily_dailyfiles/"
assar_dirs[out_daily_yearlyfiles]=${assar_dirs[hpc_data]}"mrms_nc_preciprate_daily_yearlyfiles/"
assar_dirs[out_zip_fullres_dailyfiles]=${assar_dirs[hpc_data]}"mrms_nc_preciprate_fullres_yearlyfiles_zipped/"
assar_dirs[out_yearly_singlefile]=${assar_dirs[hpc_data]}"mrms_nc_preciprate_yearly_singlefile.nc"
assar_dirs[out_yearly_singlefile_stageiv]=${assar_dirs[hpc_data]}"stageiv_nc_preciprate_yearly_singlefile.nc"
assar_dirs[out_fullres_yearly_atgages]=${assar_dirs[hpc_data]}"mrms_nc_preciprate_fullres_yearlyfiles_atgages/"

#### outputs: plots
assar_dirs[out_plots_h]=${assar_dirs[repo]}"plots/h_annual_statistics/"

#### outputs: CSVs
assar_dirs[out_fullres_yearly_csvs_atgages]=${assar_dirs[hpc_data]}"mrms_csv_preciprate_fullres_yearlyfiles_atgages/"

### input datasets
assar_dirs[shp_states]=${assar_dirs[hpc_data]}"geospatial/States_shapefile.shp"
assar_dirs[shp_gages]=${assar_dirs[hpc_data]}"geospatial/rain_gages.shp"