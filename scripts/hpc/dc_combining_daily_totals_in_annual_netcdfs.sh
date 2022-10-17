#!/bin/bash
#SBATCH -D /project/quinnlab/dcl3nd/norfolk		 # working directory
#SBATCH -o /project/quinnlab/dcl3nd/norfolk/scripts/script_out_dc/job.d1.%j.%N.out	# Name of the output file (eg. myMPI.oJobID)
#SBATCH --ntasks=1				# Number of tasks per serial job (must be 1)
#SBATCH -p standard				# Queue name "standard" (serial)
#SBATCH -A quinnlab_paid				# allocation name
#SBATCH -t 36:00:00				# Run time per serial job (hh:mm:ss)
#SBATCH --array=1-22			# Array of jobs to loop through 22 years (2001-2022)
#SBATCH --mem-per-cpu=30000

# sacct -j <jobid> -o jobid,jobname%20,user,partition,state,start,end,NodeList%60
# d4b_combining_daily_totals_in_annual_netcdfs.sh

module purge
module load anaconda
source activate mrms_processing

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
python scripts/_d4b_combining_daily_totals_in_annual_netcdfs.py ${year} "/project/quinnlab/dcl3nd/norfolk/data/mrms_for_rainyday_daily/" "/project/quinnlab/dcl3nd/norfolk/data/_scratch_zarrs/" "/project/quinnlab/dcl3nd/norfolk/data/mrms_for_rainyday_daily_consolidated/"


