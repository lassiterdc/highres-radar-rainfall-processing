#!/bin/bash
#SBATCH -o _script_outputs/%x/%A_%a_%N.out
#SBATCH -e _script_errors/%x/%A_%a_%N.out
#SBATCH --ntasks=1				# Number of tasks per serial job (must be 1)
#SBATCH -p standard				# Queue name "standard" (serial)
#SBATCH -A quinnlab				# allocation name
#SBATCH -t 24:00:00				# Run time per serial job (hh:mm:ss)
#SBATCH --array=2-3			#  1-3 (currently doing 3 groupings)
#SBATCH --mem-per-cpu=32000
#SBATCH --mail-user=dcl3nd@virginia.edu          # address for email notification
#SBATCH --mail-type=ALL   
# SBATCH --exclude=udc-ba26-18,udc-ba27-14,udc-ba26-16,udc-ba26-17

# module purge
# module load gcc openmpi eccodes anaconda
# conda create -n mrms_qaqc
# source activate mrms_qaqc
# conda install -c conda-forge xarray dask netCDF4 bottleneck h5netcdf zarr flox scipy cftime

# ijob -c 1 -A quinnlab_paid -p standard --time=0-08:00:00 --mem-per-cpu=32000

# created a new environment for this script due to a weird error

module purge
module load gcc openmpi eccodes anaconda # the stuff other than anaconda was to ensure eccodes loaded correctly
source activate mrms_qaqc

source __utils.sh
source __directories.sh

python ${assar_dirs[hpc_da3b]} ${assar_dirs[out_fullres_dailyfiles_consolidated]} ${SLURM_ARRAY_TASK_ID}


