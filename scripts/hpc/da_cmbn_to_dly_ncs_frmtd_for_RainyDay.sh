#!/bin/bash
#SBATCH -o _script_outputs/%x/%A_%a_%N.out
#SBATCH -e _script_errors/%x/%A_%a_%N.out
#SBATCH --ntasks=1				# Number of tasks per serial job (must be 1)
#SBATCH -p standard				# Queue name "standard" (serial)
#SBATCH -A quinnlab_paid				# allocation name
#SBATCH -t 48:00:00				# Run time per serial job (hh:mm:ss)
#SBATCH --array=10,11,12,3,4,5,7,8,9			#  1-366 Array of jobs to loop through (366 days)
#SBATCH --mem-per-cpu=80000
#SBATCH --mail-user=dcl3nd@virginia.edu          # address for email notification
#SBATCH --mail-type=ALL   
#SBATCH --exclude=udc-ba26-18,udc-ba27-14,udc-ba26-16,udc-ba26-17

module purge
module load gcc openmpi eccodes anaconda # the stuff other than anaconda was to ensure eccodes loaded correctly
source activate mrms_processing

source __utils.sh
source __directories.sh
# move to working directory
cd ${assar_dirs[repo]}

# all years, hours and minutes to loop through for each day of the year
YEARS=$(seq 2001 2022)

# loop through all years
for YEAR in ${YEARS}
do
	year=${YEAR}
	determine_month_and_day ${YEAR} ${SLURM_ARRAY_TASK_ID}
	month=${array_out[0]}
	day=${array_out[1]}
	# process the mrms mesonet grib data
	# echo "Node ID: $HOSTNAME"
	# echo "Slurm Array Task ID: ${SLURM_ARRAY_TASK_ID}"
	python ${assar_dirs[hpc_da]}  ${year}${month}${day} ${assar_dirs[raw_mrms]} ${assar_dirs[raw_nssl]} ${assar_dirs[raw_mrms_quantized]} ${assar_dirs[scratch_zarrs]} ${assar_dirs[scratch_gribs]} ${assar_dirs[out_fullres_dailyfiles]}
	# echo "Finished attempt to create netcdf for ${year}${month}${day}"
done

