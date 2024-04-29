#!/bin/bash
#SBATCH -o _script_outputs/%x/%A_%a_%N.out
#SBATCH -e _script_errors/%x/%A_%a_%N.out
#SBATCH --ntasks=1				# Number of tasks per serial job (must be 1)
#SBATCH -p standard				# Queue name "standard" (serial)
#SBATCH -A quinnlab				# allocation name
#SBATCH -t 10:00:00				# Run time per serial job (hh:mm:ss)
#SBATCH --mem=50000
#SBATCH --array=2-11,15-22			# Array of jobs to loop through (2002-2011, 2015-2021)
#SBATCH --mail-user=dcl3nd@virginia.edu          # address for email notification
#SBATCH --mail-type=ALL   

# ijob -c 1 -A quinnlab -p standard --time=0-08:00:00 --mem-per-cpu=36000

module purge
module load gcc openmpi eccodes anaconda
source activate mrms_analysis

source __directories.sh
# move to working directory
# cd ${assar_dirs[repo]}

if [ ${SLURM_ARRAY_TASK_ID} -lt 10 ]
then
  year=200${SLURM_ARRAY_TASK_ID}
else
  year=20${SLURM_ARRAY_TASK_ID}
fi

# echo "Running python script to create csv at gages and netcdf spanning gages for year ${year}..."
python ${assar_dirs[hpc_i]} "${assar_dirs[out_fullres_dailyfiles_consolidated]}${year}*.nc" "${assar_dirs[out_fullres_dailyfiles]}${year}*.nc" "${assar_dirs[stageiv_rainfall]}${year}/*.nc" ${assar_dirs[shp_gages]} "${assar_dirs[out_fullres_yearly_csvs_atgages]}${year}.csv" "${assar_dirs[out_fullres_yearly_atgages]}${year}.nc" "${assar_dirs[out_fullres_yearly_csvs_atgages_stageiv]}${year}.csv" "${assar_dirs[out_fullres_yearly_atgages_stageiv]}${year}.nc"
# echo "Script complete."