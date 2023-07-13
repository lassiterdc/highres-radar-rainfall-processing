#!/bin/bash
#SBATCH -o _script_outputs/%x/%A_%a_%N.out
#SBATCH -e _script_errors/%x/%A_%a_%N.out
#SBATCH --ntasks=1				# Number of tasks per serial job (must be 1)
#SBATCH -p standard				# Queue name "standard" (serial)
#SBATCH -A quinnlab_paid				# allocation name
#SBATCH -t 1:00:00				# Run time per serial job (hh:mm:ss)
#SBATCH --array=1			#  1-366 Array of jobs to loop through (366 days)
# SBATCH --mem-per-cpu=80000
#SBATCH --mail-user=dcl3nd@virginia.edu          # address for email notification
#SBATCH --mail-type=ALL   
# SBATCH --exclude=udc-ba26-18,udc-ba27-14,udc-ba26-16,udc-ba26-17

module purge
module load gcc openmpi eccodes anaconda # the stuff other than anaconda was to ensure eccodes loaded correctly
source activate mrms_processing

source __utils.sh
source __directories.sh

python ${assar_dirs[hpc_da3]} ${assar_dirs[out_fullres_dailyfiles]} ${assar_dirs[out_fullres_dailyfiles_consolidated]} ${assar_dirs[scratch_zarrs]} ${assar_dirs[scratch_csv]}


