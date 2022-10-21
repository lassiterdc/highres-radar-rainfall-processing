#!/bin/bash
#SBATCH -o _script_outputs/%x.out
#SBATCH -e _script_errors/%x.out
#SBATCH --ntasks=1				# Number of tasks per serial job (must be 1)
#SBATCH -p standard				# Queue name "standard" (serial)
#SBATCH -A quinnlab_paid				# allocation name
#SBATCH -t 36:00:00				# Run time per serial job (hh:mm:ss)
#SBATCH --array=1-22			# Array of jobs to loop through 22 years (2001-2022)
#SBATCH --mem-per-cpu=30000

module purge
module load anaconda
source activate mrms_processing

source _work__utils.sh
source __directories.sh
# move to working directory
cd ${assar_dirs[repo]}

# assign year
if [ ${SLURM_ARRAY_TASK_ID} -lt 10 ]
then
	year=200${SLURM_ARRAY_TASK_ID}
else
	year=20${SLURM_ARRAY_TASK_ID}
fi

echo "Node ID: $HOSTNAME"
echo "Slurm Array Task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Attempting to combine daily total netcdfs into single netcdf for year ${year}"
python ${assar_dirs[hpc_dc]} ${year} ${assar_dirs[out_daily_dailyfiles]} ${assar_dirs[scratch_zarrs]} ${assar_dirs[out_daily_yearlyfiles]}


