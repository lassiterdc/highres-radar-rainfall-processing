#!/bin/bash
#SBATCH -o _script_outputs/%x/%A_%a_%N.out
#SBATCH -e _script_errors/%x/%A_%a_%N.out
#SBATCH --ntasks=1				# Number of tasks per serial job (must be 1)
#SBATCH -p standard				# Queue name "standard" (serial)
#SBATCH -A quinnlab_paid				# allocation name
#SBATCH -t 10:00:00				# Run time per serial job (hh:mm:ss)
#SBATCH --mem=50000
#SBATCH --array=1
#SBATCH --mail-user=dcl3nd@virginia.edu          # address for email notification
#SBATCH --mail-type=ALL   

module purge
module load gcc openmpi eccodes anaconda
# source activate xesmf_env
source activate mrms_analysis

source __directories.sh

# echo "Running script to generate plots of annual totals and annual anomolies..."
python ${assar_dirs[hpc_hb2]} ${assar_dirs[out_yearly_singlefile]} ${assar_dirs[out_yearly_singlefile_stageiv]} ${assar_dirs[shp_states]} ${assar_dirs[out_plots_h]}
# echo "Script finished."