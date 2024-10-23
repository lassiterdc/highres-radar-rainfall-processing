#!/bin/bash
#SBATCH -o _script_outputs/%x/%A_out_%a_%N.out
#SBATCH -e _script_outputs/%x/%A_err_%a_%N.out
#SBATCH --ntasks=1				# Number of tasks per serial job (must be 1)
#SBATCH -p standard				# Queue name "standard" (serial)
#SBATCH -A quinnlab				# allocation name
#SBATCH -t 24:00:00				# Run time per serial job (hh:mm:ss)
#SBATCH --array=1-366			#  1-366 Array of jobs to loop through (366 days)
#SBATCH -c 1
# SBATCH --mem-per-cpu=36000
#SBATCH --mail-user=dcl3nd@virginia.edu          # address for email notification
#SBATCH --mail-type=ALL   

# cd /project/quinnlab/dcl3nd/norfolk/highres-radar-rainfall-processing/scripts/hpc
# ijob -A quinnlab -p standard --time=0-08:00:00 -c 1

module purge
source __directories.sh
source __utils.sh

dir_outs=_script_outputs/
mkdir -p ${dir_outs}${SLURM_JOB_NAME}
archive_previous_script_outfiles

module load gcc openmpi eccodes miniforge
DIR=~/.conda/envs/rainyday
source activate mrms_processing
export PATH=$DIR/bin:$PATH
export LD_LIBRARY_PATH=$DIR/lib:$PATH
export PYTHONPATH=$DIR/lib/python3.11/site-packages:$PATH

# SLURM_ARRAY_TASK_ID=36

YEARS=$(seq 2001 2024)

# loop through all years
for YEAR in ${YEARS}
do
	year=${YEAR}
	determine_month_and_day ${YEAR} ${SLURM_ARRAY_TASK_ID}
	month=${array_out[0]}
	day=${array_out[1]}
	rm -rf ${assar_dirs[out_fullres_dailyfiles_consolidated]}*${year}${month}${day}*
	rm -rf ${assar_dirs[scratch_zarrs]}*${year}${month}${day}*
	# other locs I was temporarily storing stuff
	rm -rf /project/quinnlab/dcl3nd/norfolk/highres-radar-rainfall-processing/data/_scratch/zarrs_delete/*${year}${month}${day}*
	rm -rf /scratch/dcl3nd/highres-radar-rainfall-processing/_scratch/zarrs_delete/*${year}${month}${day}*.zarr
	rm -rf /scratch/dcl3nd/highres-radar-rainfall-processing/_scratch/zarrs_delete_2/*${year}${month}${day}*
done