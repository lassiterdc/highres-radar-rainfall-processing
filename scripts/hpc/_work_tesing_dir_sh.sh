#!/bin/bash
source __directories.sh
FILE_OUT="${assar_dirs[working]}${assar_dirs[out_test]}"

#SBATCH -D ${assar_dirs[working]}	 # working directory
#SBATCH -o $FILE_OUT	# Name of the output file (eg. myMPI.oJobID)
#SBATCH --ntasks=1				# Number of tasks per serial job (must be 1)
#SBATCH -p standard				# Queue name "standard" (serial)
#SBATCH -A quinnlab_paid			# allocation name
#SBATCH -t 1:00:00				# Run time per serial job (hh:mm:ss)
#SBATCH --array=1-11			# Array of jobs to loop through

# sample download link
  # https://griffin-objstore.opensciencedatacloud.org/noaa-mrms-reanalysis/MRMS_PrecipRate_2001.tar

if [ ${SLURM_ARRAY_TASK_ID} -lt 10 ]
then
  year=200${SLURM_ARRAY_TASK_ID}
else
  year=20${SLURM_ARRAY_TASK_ID}
fi

echo "Tested for ${year}"
