#!/bin/bash
#SBATCH -o _script_outputs/%x/%A_%a_%N.out
#SBATCH -e _script_errors/%x/%A_%a_%N.out
#SBATCH --ntasks=1				# Number of tasks per serial job (must be 1)
#SBATCH -p standard				# Queue name "standard" (serial)
#SBATCH -A quinnlab_paid				# allocation name
#SBATCH -t 72:00:00				# Run time per serial job (hh:mm:ss)
#SBATCH --array=1-366			# Array of jobs to loop through (366 days)
#SBATCH --mail-user=dcl3nd@virginia.edu          # address for email notification
#SBATCH --mail-type=ALL   

# every x minutes until this node's task is complete:
    # return some value indicating whether the node is complete
        # if true
            # stop this script
        # else
            # continue
    # call python script to compute median runtime for each task echo'd in the .out file
    # call python script to compute median runtime for each task echo'd in ALL .out files
    # if median runtime for THIS node is significantly more than median runtime for all nodes
        # reset this job on a different node
