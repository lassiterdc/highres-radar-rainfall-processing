#!/bin/bash
#SBATCH -o _script_outputs/%x.out
#SBATCH -e _script_errors/%x.out
#SBATCH --ntasks=1				# Number of tasks per serial job (must be 1)
#SBATCH -p standard				# Queue name "standard" (serial)
#SBATCH -A quinnlab_paid			# allocation name
#SBATCH -t 10:00:00				# Run time per serial job (hh:mm:ss)
#SBATCH --array=1-11			# Array of jobs to loop through

# sample download link
  # https://griffin-objstore.opensciencedatacloud.org/noaa-mrms-reanalysis/MRMS_PrecipRate_2001.tar

source __utils.sh
source __directories.sh
# confirm working directory exists
mkdir -p ${assar_dirs[repo]}${assar_dirs[raw_nssl]}
# move to working directory
cd ${assar_dirs[repo]}${assar_dirs[raw_nssl]}

if [ ${SLURM_ARRAY_TASK_ID} -lt 10 ]
then
  year=200${SLURM_ARRAY_TASK_ID}
else
  year=20${SLURM_ARRAY_TASK_ID}
fi

wget -q -c https://griffin-objstore.opensciencedatacloud.org/noaa-mrms-reanalysis/MRMS_PrecipRate_${year}.tar
# echo "downloaded data for year $year"
# https://www.pendrivelinux.com/how-to-open-a-tar-file-in-unix-or-linux/
tar -xf MRMS_PrecipRate_${year}.tar
# echo "unzipped tar file for year $year"
rm MRMS_PrecipRate_${year}.tar

echo "Finished downloading and unzipping data for year ${year}"
