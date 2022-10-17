#!/bin/bash
#SBATCH -D /gpfs/gpfs0/project/quinnlab/dcl3nd/norfolk/scripts/		 # working directory
#SBATCH -o /project/quinnlab/dcl3nd/norfolk/scripts/script_out_e/job.nssl_at_gages.%j.%N.out	# Name of the output file (eg. myMPI.oJobID)
#SBATCH --ntasks=1				# Number of tasks per serial job (must be 1)
#SBATCH -p standard				# Queue name "standard" (serial)
#SBATCH -A quinnlab_paid				# allocation name
#SBATCH -t 10:00:00				# Run time per serial job (hh:mm:ss)
#SBATCH --array=1			# Array of jobs to loop through (366 days)

# module purge
# module load anaconda
# source activate norfolk

#!/bin/bash
# NOTE : combine output files into single file
cat /project/quinnlab/dcl3nd/norfolk/scripts/script_out_d1/*.out > "/project/quinnlab/dcl3nd/norfolk/scripts/_out_d1_cmbnd.txt"
cat /project/quinnlab/dcl3nd/norfolk/scripts/script_out_d4a/*.out > "/project/quinnlab/dcl3nd/norfolk/scripts/_out_d4a_cmbnd.txt"
cat /project/quinnlab/dcl3nd/norfolk/scripts/script_out_d4b/*.out > "/project/quinnlab/dcl3nd/norfolk/scripts/_out_d4b_cmbnd.txt"
cat /project/quinnlab/dcl3nd/norfolk/scripts/script_out_i/*.out > "/project/quinnlab/dcl3nd/norfolk/scripts/_out_i_cmbnd.txt"

# for i in {1..366}
# do
#    seff 42310329_"$i" >> "/project/quinnlab/dcl3nd/norfolk/scripts/d_seff_cmbnd.txt"
#    echo "$i"
# done