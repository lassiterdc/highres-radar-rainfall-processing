#!/bin/bash
#SBATCH -D /project/quinnlab/dcl3nd/norfolk/		 # working directory
#SBATCH -o /project/quinnlab/dcl3nd/norfolk/scripts/script_out_h/job.nssl_at_gages.%j.%N.out	# Name of the output file (eg. myMPI.oJobID)
#SBATCH --ntasks=1				# Number of tasks per serial job (must be 1)
#SBATCH -p standard				# Queue name "standard" (serial)
#SBATCH -A quinnlab_paid				# allocation name
#SBATCH -t 10:00:00				# Run time per serial job (hh:mm:ss)
#SBATCH --mem=50000
#SBATCH --array=1

module purge
module load gcc openmpi eccodes anaconda
source activate mrms_analysis

echo "Running python script to create csv at gages and netcdf spanning gages for year ${year}..."
python "scripts/_work_h_generate_annual_statistics_plots_to_qaqc.py"  "data/mrms_for_rainyday/${year}*.nc" "data/shapefiles/rain_gages.shp" "data/mrms_for_rainyday_subset_norfolk_csvs/${year}.csv" "data/mrms_for_rainyday_subset_norfolk_netcdfs/${year}.nc"
echo "Script complete."