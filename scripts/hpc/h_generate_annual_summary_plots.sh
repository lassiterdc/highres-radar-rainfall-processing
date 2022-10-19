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

echo "Running python script to create mrms_daily_totals.nc and mrms_yearly_totals.nc..."
python "scripts/_ha_generate_annual_statistics_netcdfs.py" "data/mrms_for_rainyday_daily_consolidated/*.nc" 'data/mrms_daily_totals.nc' 'data/mrms_yearly_totals.nc' '_scratch_zarrs/h_yearly.zarr' "data/shapefiles/States_shapefile.shp" "plots/h_annual_statistics/"
echo "First script complete. Running script to generate plots of annual totals and annual anomolies."
python "scripts/_hb_generate_annual_statistics_plots_to_qaqc.py" 'data/mrms_daily_totals.nc' 'data/mrms_yearly_totals.nc' "data/shapefiles/States_shapefile.shp" "plots/h_annual_statistics/"
echo "Script finished."