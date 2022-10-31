#!/bin/bash
#SBATCH -o _script_outputs/%x_%A.out
#SBATCH -e _script_errors/%x_%A.out
#SBATCH --ntasks=1				# Number of tasks per serial job (must be 1)
#SBATCH -p standard				# Queue name "standard" (serial)
#SBATCH -A quinnlab_paid				# allocation name
#SBATCH -t 10:00:00				# Run time per serial job (hh:mm:ss)
#SBATCH --mem=50000
#SBATCH --array=1-22			# Array of jobs to loop through (2001-2022)

module purge
module load gcc openmpi eccodes anaconda
source activate mrms_analysis

source __directories.sh
# move to working directory
cd ${assar_dirs[repo]}

if [ ${SLURM_ARRAY_TASK_ID} -lt 10 ]
then
  year=200${SLURM_ARRAY_TASK_ID}
else
  year=20${SLURM_ARRAY_TASK_ID}
fi

# echo "Running python script to create csv at gages and netcdf spanning gages for year ${year}..."
python ${assar_dirs[hpc_i]}  "${assar_dirs[out_fullres_dailyfiles]}${year}*.nc" ${assar_dirs[shp_gages]} "${assar_dirs[out_fullres_yearly_csvs_atgages]}${year}.csv" "${assar_dirs[out_fullres_yearly_atgages]}${year}.nc"
# echo "Script complete."