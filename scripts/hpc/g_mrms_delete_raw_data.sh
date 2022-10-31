#!/bin/bash
#SBATCH -o _script_outputs/%x_%A.out
#SBATCH -e _script_errors/%x_%A.out
#SBATCH --ntasks=1				# Number of tasks per serial job (must be 1)
#SBATCH -p standard				# Queue name "standard" (serial)
#SBATCH -A quinnlab				# allocation name
#SBATCH -t 48:00:00				# Run time per serial job (hh:mm:ss)
#SBATCH --array=1-3			# Array of jobs to loop through

source __directories.sh
# move to working directory
cd ${assar_dirs[repo]}

if [ ${SLURM_ARRAY_TASK_ID} = 1 ]
then
  # echo "removing nssl gribs"
  rm -rf ${assar_dirs[raw_nssl]}
  echo "removed nssl gribs"
fi

if [ ${SLURM_ARRAY_TASK_ID} = 2 ]
then
  # echo "removing mesonet gribs"
  rm -rf ${assar_dirs[raw_mrms]}
  echo "removed mesonet gribs"
fi

if [ ${SLURM_ARRAY_TASK_ID} = 3 ]
then
  # echo "removing mesonet netcdfs"
  rm -rf ${assar_dirs[raw_mrms_quantized]}
  echo "removed mesonet netcdfs"
fi