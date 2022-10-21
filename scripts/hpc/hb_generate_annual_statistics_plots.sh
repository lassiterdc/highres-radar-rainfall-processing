#!/bin/bash
#SBATCH -o _script_outputs/%x.out
#SBATCH -e _script_errors/%x.out
#SBATCH --ntasks=1				# Number of tasks per serial job (must be 1)
#SBATCH -p standard				# Queue name "standard" (serial)
#SBATCH -A quinnlab_paid				# allocation name
#SBATCH -t 10:00:00				# Run time per serial job (hh:mm:ss)
#SBATCH --mem=50000
#SBATCH --array=1

module purge
module load gcc openmpi eccodes anaconda
source activate mrms_analysis

source __directories.sh
# move to working directory
cd ${assar_dirs[repo]}

echo "Running script to generate plots of annual totals and annual anomolies..."
python ${assar_dirs[hpc_hb]} ${assar_dirs[out_yearly_singlefile]} ${assar_dirs[shp_states]} ${assar_dirs[out_plots_h]}
echo "Script finished."