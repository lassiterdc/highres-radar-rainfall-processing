#!/bin/bash
#SBATCH -o _script_outputs/%x/%A_%a_%N.out
#SBATCH -e _script_errors/%x/%A_%a_%N.out
#SBATCH --ntasks=1				# Number of tasks per serial job (must be 1)
#SBATCH -p standard				# Queue name "standard" (serial)
#SBATCH -A quinnlab				# allocation name
#SBATCH -t 72:00:00				# Run time per serial job (hh:mm:ss)
#SBATCH --array=1-366			#  1-366 Array of jobs to loop through (366 days)
#SBATCH -c 1
#SBATCH --mem-per-cpu=9000
#SBATCH --mail-user=dcl3nd@virginia.edu          # address for email notification
#SBATCH --mail-type=ALL   
# SBATCH --exclude=udc-ba26-18,udc-ba27-14,udc-ba26-16,udc-ba26-17

# ijob -A quinnlab -p standard --time=0-08:00:00 -c 1 --mem-per-cpu=9000

# ijob -A quinnlab -p standard --time=0-08:00:00 -c 1 --mem-per-cpu=80000 | ymd=20200406 | 8.51 minutes
# ijob -A quinnlab -p standard --time=0-08:00:00 -c 4 --mem-per-cpu=20000 | ymd=20210406 | 8.87 minutes
# ijob -A quinnlab -p standard --time=0-08:00:00 -c 8 --mem-per-cpu=9000 | ymd=20220406 | 
# ijob -A quinnlab -p standard --time=0-08:00:00 -c 1 --mem-per-cpu=9000 | ymd=20230406 | 20 minutes
# python ${assar_dirs[hpc_da]} $ymd ${assar_dirs[raw_mrms]} ${assar_dirs[raw_nssl]} ${assar_dirs[raw_mrms_quantized]} ${assar_dirs[scratch_zarrs]} ${assar_dirs[scratch_gribs]} ${assar_dirs[out_fullres_dailyfiles]}

# -c 1 --mem-per-cpu=9000 | 15 minutes
# --time=0-08:00:00 -c 1 --mem-per-cpu=80000, SLURM_ARRAY_TASK_ID=105, YEAR = 2023 exporting to zarr | 14.91 minutes
# made some changes (Exporting zarr earlier in process)
# -c 2 --mem-per-cpu=80000 took forever


# inspecting errors
# cd /project/quinnlab/dcl3nd/norfolk/highres-radar-rainfall-processing/scripts/hpc/_script_errors/da_cmbn_to_dly_ncs_frmtd_for_RainyDay.sh/
# inspecting outputs
# cd /project/quinnlab/dcl3nd/norfolk/highres-radar-rainfall-processing/scripts/hpc/_script_outputs/da_cmbn_to_dly_ncs_frmtd_for_RainyDay.sh/

# doing git pull
# rm /project/quinnlab/dcl3nd/norfolk/highres-radar-rainfall-processing/scripts/hpc/da_cmbn_to_dly_ncs_frmtd_for_RainyDay.sh
# rm /project/quinnlab/dcl3nd/norfolk/highres-radar-rainfall-processing/scripts/hpc/_da_cmbn_to_dly_ncs_frmtd_for_RainyDay.py
# git pull

# clearing outputs
# rm /project/quinnlab/dcl3nd/norfolk/highres-radar-rainfall-processing/scripts/hpc/_script_outputs/da_cmbn_to_dly_ncs_frmtd_for_RainyDay.sh/*
# rm /project/quinnlab/dcl3nd/norfolk/highres-radar-rainfall-processing/scripts/hpc/_script_errors/da_cmbn_to_dly_ncs_frmtd_for_RainyDay.sh/*

# cd /project/quinnlab/dcl3nd/norfolk/highres-radar-rainfall-processing/scripts/hpc
module purge
module load gcc openmpi eccodes anaconda # the stuff other than anaconda was to ensure eccodes loaded correctly
source activate mrms_processing

source __utils.sh
source __directories.sh

YEARS=$(seq 2001 2024)
# loop through all years
for YEAR in ${YEARS}
do
	determine_month_and_day ${YEAR} ${SLURM_ARRAY_TASK_ID}
	month=${array_out[0]}
	day=${array_out[1]}
	# process the mrms mesonet grib data
	# echo "Node ID: $HOSTNAME"
	# echo "Slurm Array Task ID: ${SLURM_ARRAY_TASK_ID}"
	python ${assar_dirs[hpc_da]} ${YEAR}${month}${day} ${assar_dirs[raw_mrms]} ${assar_dirs[raw_nssl]} ${assar_dirs[raw_mrms_quantized]} ${assar_dirs[scratch_zarrs]} ${assar_dirs[scratch_gribs]} ${assar_dirs[out_fullres_dailyfiles]}
	# echo "Finished attempt to create netcdf for ${YEAR}${month}${day}"
done

