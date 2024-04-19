#!/bin/bash
#SBATCH -o _script_outputs/%x/%A_%a_%N.out
#SBATCH -e _script_errors/%x/%A_%a_%N.out
#SBATCH --ntasks=1				# Number of tasks per serial job (must be 1)
#SBATCH -p standard				# Queue name "standard" (serial)
#SBATCH -A quinnlab				# allocation name
#SBATCH -t 36:00:00				# Run time per serial job (hh:mm:ss)
#SBATCH --array=1-366			# Array of jobs to loop through (366 days)
#SBATCH --mem-per-cpu=60000
#SBATCH --mail-user=dcl3nd@virginia.edu          # address for email notification
#SBATCH --mail-type=ALL   

# ijob -c 1 -A quinnlab -p standard --time=0-06:00:00 --mem-per-cpu=30000

source __utils.sh
source __directories.sh
module purge
module load anaconda
DIR=~/.conda/envs/rainyday
source activate mrms_processing
export PATH=$DIR/bin:$PATH
export LD_LIBRARY_PATH=$DIR/lib:$PATH
export PYTHONPATH=$DIR/lib/python3.11/site-packages:$PATH

# move to working directory
# cd ${assar_dirs[repo]}

# all years, hours and minutes to loop through for each day of the year
YEARS=$(seq 2002 2021)

# loop through all years
for YEAR in ${YEARS}
do
	year=${YEAR}
	determine_month_and_day ${YEAR} ${SLURM_ARRAY_TASK_ID}
	month=${array_out[0]}
	day=${array_out[1]}
#  echo "Node ID: $HOSTNAME"
#  echo "Slurm Array Task ID: ${SLURM_ARRAY_TASK_ID}"
#  echo "Attempting to create hourly and daily netcdfs for ${year}${month}${day}"
  python ${assar_dirs[hpc_db]} ${year}${month}${day} ${assar_dirs[out_fullres_dailyfiles_consolidated]} ${assar_dirs[scratch_zarrs]} ${assar_dirs[out_hourly_dailyfiles]} ${assar_dirs[out_daily_dailyfiles]} ${assar_dirs[scratch_csv]}
done