#!/bin/bash
#SBATCH -o _script_outputs/%x/%A_out_%a_%N.out
#SBATCH -e _script_outputs/%x/%A_err_%a_%N.out
#SBATCH --ntasks=1
#SBATCH -p standard
#SBATCH -A quinnlab
#SBATCH -t 5:30:00
#SBATCH --array=1-24
#SBATCH --mail-user=dcl3nd@virginia.edu
#SBATCH --mail-type=ALL 

# metadata: https://registry.opendata.aws/noaa-nws-aorc/

mkdir -p _slurm_outputs/${SLURM_JOB_NAME}
source __utils.sh
archive_previous_script_outfiles

# ijob -c 1 -A quinnlab -p standard --time=0-08:00:00
# cd /project/quinnlab/dcl3nd/norfolk/highres-radar-rainfall-processing/scripts/hpc
# SLURM_ARRAY_TASK_ID=23

module purge
source __directories.sh
module load awscli
# make sure directory exists
mkdir -p ${assar_dirs[raw_aorc]}

# move into directory
cd ${assar_dirs[raw_aorc]}

START_YEAR=2001
YEAR=$(($START_YEAR + $SLURM_ARRAY_TASK_ID - 1))

# Print the year for debugging
echo "Downloading data for year: $YEAR"
START_TIME=$(date +%s)  # Record start time

# Calculate the total size of the Zarr directory (in bytes)
TOTAL_SIZE=$(aws s3 ls --no-sign-request s3://noaa-nws-aorc-v1-1-1km/${YEAR}.zarr/ --recursive | awk '{sum += $3} END {print sum}')
HUMAN_SIZE=$(numfmt --to=iec-i --suffix=B $TOTAL_SIZE)

echo "Total size for ${YEAR}.zarr: $HUMAN_SIZE"

# Base URL for the AORC dataset in Zarr format
while true; do
    du -sh ./data/${YEAR}.zarr/ 2>/dev/null || echo "Waiting for data to be downloaded..."
    sleep 60  # Print size every 60 seconds
done &
MONITOR_PID=$!
aws s3 sync --no-sign-request s3://noaa-nws-aorc-v1-1-1km/${YEAR}.zarr/ ./data/${YEAR}.zarr/ --quiet

# stop monitoring
kill $MONITOR_PID 2>/dev/null

# Calculate total elapsed time
END_TIME=$(date +%s)
ELAPSED_TIME=$(($END_TIME - $START_TIME))

# Print the total download time in a readable format (hours:minutes:seconds)
printf "Download for year %d completed in %02d:%02d:%02d\n" \
  "$YEAR" $(($ELAPSED_TIME / 3600)) $(( ($ELAPSED_TIME % 3600) / 60 )) $(($ELAPSED_TIME % 60))

