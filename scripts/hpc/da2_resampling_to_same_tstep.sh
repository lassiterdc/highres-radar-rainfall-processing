#!/bin/bash
#SBATCH -o _script_outputs/%x/%A_output_%a_%N.out
#SBATCH -e _script_outputs/%x/%A_error_%a_%N.out
#SBATCH --ntasks=1				# Number of tasks per serial job (must be 1)
#SBATCH -p standard				# Queue name "standard" (serial)
#SBATCH -A quinnlab				# allocation name
#SBATCH -t 24:00:00				# Run time per serial job (hh:mm:ss)
#SBATCH --array=1-366			#  1-366 Array of jobs to loop through (366 days)
#SBATCH -c 4
#SBATCH --mem-per-cpu=9000
#SBATCH --mail-user=dcl3nd@virginia.edu          # address for email notification
#SBATCH --mail-type=ALL   

# SBATCH --exclude=udc-ba26-18,udc-ba27-14,udc-ba26-16,udc-ba26-17

mkdir -p -p _script_errors/${SLURM_JOB_NAME}
# ijob -A quinnlab -p standard --time=0-08:00:00 -c 4 --mem-per-cpu=9000 # | yearmonthday=20210724
# ijob -A quinnlab -p standard --time=0-08:00:00 -c 1 --mem-per-cpu=9000 # | yearmonthday=20210723
# ijob -A quinnlab -p standard --time=0-08:00:00 -c 1 --mem-per-cpu=80000 # | yearmonthday=20210722
# ijob -A quinnlab -p standard --time=0-08:00:00 -c 4 --mem-per-cpu=18000 # | yearmonthday=20210721
# ijob -A quinnlab -p standard --time=0-08:00:00 -c 1 --mem-per-cpu=36000 # | yearmonthday=20210720


# year=2021
# month=07
# day=24

# cd /project/quinnlab/dcl3nd/norfolk/highres-radar-rainfall-processing/scripts/hpc

# rm /project/quinnlab/dcl3nd/norfolk/highres-radar-rainfall-processing/scripts/hpc/_da2_resampling_to_same_tstep.py

# git pull

source __directories.sh
source __utils.sh
module purge
module load gcc openmpi eccodes anaconda
DIR=~/.conda/envs/rainyday
source activate mrms_processing
export PATH=$DIR/bin:$PATH
export LD_LIBRARY_PATH=$DIR/lib:$PATH
export PYTHONPATH=$DIR/lib/python3.11/site-packages:$PATH

# testing:
# bias correcting entire overlapping domain of stage iv and mrms
# python ${assar_dirs[hpc_da2]} ${year}${month}${day} ${assar_dirs[out_fullres_dailyfiles]} ${assar_dirs[out_fullres_dailyfiles_consolidated]} ${assar_dirs[scratch_zarrs]} ${assar_dirs[scratch_csv]} ${assar_dirs[stageiv_rainfall]}

# move to working directory
# cd ${assar_dirs[repo]}

# all years, hours and minutes to loop through for each day of the year
YEARS=$(seq 2001 2024)

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
	# bias correcting over just transposition domain
	python ${assar_dirs[hpc_da2]} ${year}${month}${day} ${assar_dirs[out_fullres_dailyfiles]} ${assar_dirs[out_fullres_dailyfiles_consolidated]} ${assar_dirs[scratch_zarrs]} ${assar_dirs[scratch_csv]} ${assar_dirs[stageiv_rainfall]} ${assar_dirs[shp_transposition_domain]}
	# echo "Finished attempt to create netcdf for ${year}${month}${day}"
done