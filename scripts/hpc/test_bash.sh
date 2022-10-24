#!/bin/bash
#SBATCH -o _script_outputs/_archive/%x_%A.out
#SBATCH -e _script_errors/_archive/%x_%A.out
#SBATCH --ntasks=1				# Number of tasks per serial job (must be 1)
#SBATCH -p standard				# Queue name "standard" (serial)
#SBATCH -A quinnlab_paid			# allocation name
#SBATCH -t 00:05:00				# Run time per serial job (hh:mm:ss)
#SBATCH --array=1-11			# Array of jobs to loop through


if [ ${SLURM_ARRAY_TASK_ID} -lt 10 ]
then
  year=200${SLURM_ARRAY_TASK_ID}
else
  year=20${SLURM_ARRAY_TASK_ID}
fi

echo "Finished downloading and unzipping data for year ${year}"
