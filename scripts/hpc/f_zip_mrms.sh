#!/bin/bash
#SBATCH -o _script_outputs/%x/%A_%j_%N.out
#SBATCH -e _script_errors/%x/%A_%j_%N.out
#SBATCH --ntasks=1				# Number of tasks per serial job (must be 1)
#SBATCH -p standard				# Queue name "standard" (serial)
#SBATCH -A quinnlab				# allocation name
#SBATCH -t 24:00:00				# Run time per serial job (hh:mm:ss)
#SBATCH --array=1-21			# Array of jobs to loop through
#SBATCH --mail-user=dcl3nd@virginia.edu          # address for email notification
#SBATCH --mail-type=ALL   

# sample download link
  # https://griffin-objstore.opensciencedatacloud.org/noaa-mrms-reanalysis/MRMS_PrecipRate_2001.tar

source __directories.sh
# move to working directory
cd ${assar_dirs[repo]}

if [ ${SLURM_ARRAY_TASK_ID} -lt 10 ]
then
  year=200${SLURM_ARRAY_TASK_ID}
else
  year=20${SLURM_ARRAY_TASK_ID}
fi

# echo "zipping ${year} data..."

tar cvfz "${assar_dirs[out_zip_fullres_dailyfiles]}/${year}.tar.gz" "${assar_dirs[out_fullres_dailyfiles]}/${year}*"

echo "zipped ${year} data"