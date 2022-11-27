#!/bin/bash

#%% code examples
# example to call array
# echo ${assar_dirs[working]}
# concatenate
#  echo "${assar_dirs[working]}${assar_dirs[raw_nssl]}"

# declare associative array
declare -A assar_dirs

# populate associative array with folder and filepaths; this should be the only thing that has to be changed
assar_dirs[repo]="/project/quinnlab/dcl3nd/norfolk/highres-radar-rainfall-processing/"

## paths relative to repo
### hpc script output folder
assar_dirs[hpc_outputs]="scripts/hpc/_script_outputs/"
assar_dirs[hpc_errors]="scripts/hpc/_script_errors/"

### raw data
assar_dirs[raw_nssl]="data/raw_data/mrms_grib_nssl/"
assar_dirs[raw_mrms_quantized]="data/raw_data/mrms_nc_quant/"
assar_dirs[raw_mrms]="data/raw_data/mrms_grib_mesonet/"

### scripts
assar_dirs[hpc_da]="scripts/hpc/_da_cmbn_to_dly_ncs_frmtd_for_RainyDay.py"
assar_dirs[hpc_db]="scripts/hpc/_db_resampling_to_hourly_and_daily_timesteps.py"
assar_dirs[hpc_dc]="scripts/hpc/_dc_combining_daily_totals_in_annual_netcdfs.py"
assar_dirs[hpc_ha]="scripts/hpc/_ha_generate_annual_statistics_netcdfs.py"
assar_dirs[hpc_hb]="scripts/hpc/_hb_generate_annual_statistics_plots.py"
assar_dirs[hpc_i]="scripts/hpc/_i_extract_mrms_at_gages.py"

### scratch folders
assar_dirs[scratch_zarrs]="data/_scratch/zarrs/"
assar_dirs[scratch_gribs]="data/_scratch/gribs/"

### Outputs:
#### outputs: processed data (format is [source]_[type]_[data]_[timestep]_[file aggregation]_[other qualifiers])
assar_dirs[out_fullres_dailyfiles]="data/mrms_nc_preciprate_fullres_dailyfiles/"
assar_dirs[out_hourly_dailyfiles]="data/mrms_nc_preciprate_hourly_dailyfiles/"
assar_dirs[out_daily_dailyfiles]="data/mrms_nc_preciprate_daily_dailyfiles/"
assar_dirs[out_daily_yearlyfiles]="data/mrms_nc_preciprate_daily_yearlyfiles/"
assar_dirs[out_zip_fullres_dailyfiles]="data/mrms_nc_preciprate_fullres_yearlyfiles_zipped/"
assar_dirs[out_yearly_singlefile]="data/mrms_nc_preciprate_yearly_singlefile.nc"
assar_dirs[out_fullres_yearly_atgages]="data/mrms_nc_preciprate_fullres_yearlyfiles_subset-over-gage-network/"

#### outputs: plots
assar_dirs[out_plots_h]="plots/h_annual_statistics/"

#### outputs: CSVs
assar_dirs[out_fullres_yearly_csvs_atgages]="data/mrms_csv_preciprate_fullres_yearlyfiles_subset-over-gage-network/"

### input datasets
assar_dirs[shp_states]="data/geospatial/States_shapefile.shp"
assar_dirs[shp_gages]="data/geospatial/rain_gages.shp"