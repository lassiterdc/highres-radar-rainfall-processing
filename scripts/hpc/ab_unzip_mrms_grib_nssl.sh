#!/bin/bash
#SBATCH -o _script_outputs/%x/%A_%a_%N.out
#SBATCH -e _script_errors/%x/%A_%a_%N.out
#SBATCH --ntasks=1				# Number of tasks per serial job (must be 1)
#SBATCH -p standard				# Queue name "standard" (serial)
#SBATCH -A quinnlab				# allocation name
#SBATCH -t 72:00:00				# Run time per serial job (hh:mm:ss)
#SBATCH --array=1-366			# Array of jobs to loop through (366 days)
#SBATCH --mail-user=dcl3nd@virginia.edu          # address for email notification
#SBATCH --mail-type=ALL   

# ijob -c 1 -A quinnlab -p standard --time=0-08:00:00

# SLURM_ARRAY_TASK_ID=171
# year=2011

# cd /project/quinnlab/dcl3nd/norfolk/highres-radar-rainfall-processing/scripts/hpc

# if doing a git pull, might have to delete the file before pulling
# rm /project/quinnlab/dcl3nd/norfolk/highres-radar-rainfall-processing/scripts/hpc/ab_unzip_mrms_grib_nssl.sh
# git pull

source __utils.sh
source __directories.sh
#confirm working directory exists
mkdir -p ${assar_dirs[raw_nssl]}
# move to working directory
cd ${assar_dirs[raw_nssl]}

# all years, hours and minutes to loop through for each day of the year
YEARS=$(seq 2001 2011)
for YEAR in ${YEARS}; do
SECONDS=0
year=${YEAR}
# define day's worth of .gz files to unzip
determine_month_and_day ${YEAR} ${SLURM_ARRAY_TASK_ID}
month=${array_out[0]}
day=${array_out[1]}
DATE=${year}${month}${day}
FPATTERN=*"${DATE}"*".gz"
# Create a list of .gz files to unzip given file pattern
FILES=$(ls ${FPATTERN} 2>/dev/null)
# Loop through the .gz files and unzip them
if [ -n "$FILES" ]; then
    for FILE in $FILES; do
        gunzip "$FILE"
        # echo "Unzipped ${FILE}"
    done
    duration=$SECONDS
    echo "Unzipped all .gz files for date ${month}/${day}/${year}; Time elapsed: $(($duration / 60)) minutes and $(($duration % 60)) seconds"
fi
done