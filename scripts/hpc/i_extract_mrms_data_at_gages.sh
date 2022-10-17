#!/bin/bash
#SBATCH -D /project/quinnlab/dcl3nd/norfolk/		 # working directory
#SBATCH -o /project/quinnlab/dcl3nd/norfolk/scripts/script_out_i/job.nssl_at_gages.%j.%N.out	# Name of the output file (eg. myMPI.oJobID)
#SBATCH --ntasks=1				# Number of tasks per serial job (must be 1)
#SBATCH -p standard				# Queue name "standard" (serial)
#SBATCH -A quinnlab_paid				# allocation name
#SBATCH -t 10:00:00				# Run time per serial job (hh:mm:ss)
#SBATCH --mem=50000
#SBATCH --array=1-22			# Array of jobs to loop through (2001-2022)

module purge
module load gcc openmpi eccodes anaconda
source activate mrms_analysis

if [ ${SLURM_ARRAY_TASK_ID} -lt 10 ]
then
  year=200${SLURM_ARRAY_TASK_ID}
else
  year=20${SLURM_ARRAY_TASK_ID}
fi

# FILE="data/mrms_for_rainyday_subset_norfolk_csvs/${year}.csv"

#EXISTS=0
#if [ -f "$FILE" ]; then
#  EXISTS=1
#  echo "$FILE already exists!"
#  echo "$EXISTS is the exists variable. If 1, the script should end here."
#else
#  echo "$FILE does not yet exist. Running python script to create..."
  # python "scripts/_i_extract_mrms_at_gages.py"  "data/mrms_for_rainyday/${year}*.nc" "data/shapefiles/rain_gages.shp" "data/mrms_for_rainyday_subset_norfolk_csvs/${year}.csv" "data/mrms_for_rainyday_subset_norfolk_netcdfs/${year}.nc" # "data/_scratch_zarrs/${year}_norfolk.zarr"
  # echo "$FILE created"
# fi

echo "Running python script to create csv at gages and netcdf spanning gages for year ${year}..."
python "scripts/_i_extract_mrms_at_gages.py"  "data/mrms_for_rainyday/${year}*.nc" "data/shapefiles/rain_gages.shp" "data/mrms_for_rainyday_subset_norfolk_csvs/${year}.csv" "data/mrms_for_rainyday_subset_norfolk_netcdfs/${year}.nc"
echo "Script complete."