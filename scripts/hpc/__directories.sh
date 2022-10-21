#!/bin/bash

#%% code examples
# example to call array
# echo ${assar_dirs[working]}
# concatenate
#  echo "${assar_dirs[working]}${assar_dirs[raw_nssl]}"

# declare associative array
declare -A assar_dirs

# populate associative array with folder and filepaths
assar_dirs[repo]="/project/quinnlab/dcl3nd/norfolk/highres-radar-rainfall-processing/"
## paths relative to repo
assar_dirs[raw_nssl]="data/raw_data/mrms_grib_nssl/"
assar_dirs[raw_mrms_quantized]="data/raw_data/mrms_nc_quant/"
assar_dirs[raw_mrms]="data/raw_data/mrms_grib_mesonet/"
# scripts
assar_dirs[hpc_da]="scripts/hpc/_da_cmbn_to_dly_ncs_frmtd_for_RainyDay.py"
# output folders
assar_dirs[scratch_zarrs]="data/_archive/zarrs"
assar_dirs[scratch_gribs]="data/_archive/gribs"
assar_dirs[]="data/mrms_preciprate_fullres_daily"
assar_dirs[]=""
assar_dirs[]=""
assar_dirs[]=""
assar_dirs[]=""
assar_dirs[]=""
assar_dirs[]=""
assar_dirs[]=""
assar_dirs[]=""
assar_dirs[]=""
assar_dirs[]=""
assar_dirs[]=""
assar_dirs[]=""
assar_dirs[]=""
assar_dirs[]=""
