#!/bin/bash
#SBATCH -D /scratch/dcl3nd/norfolk/data/	 # working directory
#SBATCH -o /project/quinnlab/dcl3nd/norfolk/scripts/script_out_f/job.%j.%N.out	# Name of the output file (eg. myMPI.oJobID)
#SBATCH --ntasks=1				# Number of tasks per serial job (must be 1)
#SBATCH -p standard				# Queue name "standard" (serial)
#SBATCH -A quinnlab				# allocation name
#SBATCH -t 24:00:00				# Run time per serial job (hh:mm:ss)
#SBATCH --array=1-21			# Array of jobs to loop through

# sample download link
  # https://griffin-objstore.opensciencedatacloud.org/noaa-mrms-reanalysis/MRMS_PrecipRate_2001.tar

if [ ${SLURM_ARRAY_TASK_ID} -lt 10 ]
then
  year=200${SLURM_ARRAY_TASK_ID}
else
  year=20${SLURM_ARRAY_TASK_ID}
fi

echo "zipping ${year} data..."

tar cvfz mrms_zipped/${year}.tar.gz mrms_for_rainyday/${year}*

echo "zipped ${year} data"