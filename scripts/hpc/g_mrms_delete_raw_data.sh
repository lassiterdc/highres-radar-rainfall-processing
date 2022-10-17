#!/bin/bash
#SBATCH -D /project/quinnlab/dcl3nd/norfolk/data/	 # working directory
#SBATCH -o /project/quinnlab/dcl3nd/norfolk/scripts/script_out_g/job.%j.%N.out	# Name of the output file (eg. myMPI.oJobID)
#SBATCH --ntasks=1				# Number of tasks per serial job (must be 1)
#SBATCH -p standard				# Queue name "standard" (serial)
#SBATCH -A quinnlab				# allocation name
#SBATCH -t 48:00:00				# Run time per serial job (hh:mm:ss)
#SBATCH --array=1-4			# Array of jobs to loop through

if [ ${SLURM_ARRAY_TASK_ID} = 1 ]
then
  echo "removing nssl gribs"
  rm -rf mrms_grib_nssl
  echo "removed nssl gribs"
fi

if [ ${SLURM_ARRAY_TASK_ID} = 2 ]
then
  echo "removing mesonet gribs"
  rm -rf mrms_grib_mesonet
  echo "removed mesonet gribs"
fi

if [ ${SLURM_ARRAY_TASK_ID} = 3 ]
then
  echo "removing mesonet netcdfs"
  rm -rf mrms_ncdf_mesonet
  echo "removed mesonet netcdfs"
fi

if [ ${SLURM_ARRAY_TASK_ID} = 4 ]
then
  echo "removing mesonet netcdfs from pngs"
  rm -rf mrms_netcdf_mesonet2
  echo "removed mesonet netcdfs from pngs"
fi