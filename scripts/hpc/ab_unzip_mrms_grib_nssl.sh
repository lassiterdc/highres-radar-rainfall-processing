#!/bin/bash
#SBATCH -o _script_outputs/%x/%A_%j_%N.out
#SBATCH -e _script_errors/%x/%A_%j_%N.out
#SBATCH --ntasks=1				# Number of tasks per serial job (must be 1)
#SBATCH -p standard				# Queue name "standard" (serial)
#SBATCH -A quinnlab_paid				# allocation name
#SBATCH -t 72:00:00				# Run time per serial job (hh:mm:ss)
#SBATCH --array=1-366			# Array of jobs to loop through (366 days)
#SBATCH --mail-user=dcl3nd@virginia.edu          # address for email notification
#SBATCH --mail-type=ALL   

source __utils.sh
source __directories.sh
#confirm working directory exists
mkdir -p ${assar_dirs[repo]}${assar_dirs[raw_nssl]}
# move to working directory
cd ${assar_dirs[repo]}${assar_dirs[raw_nssl]}

# all years, hours and minutes to loop through for each day of the year
YEARS=$(seq 2001 2011)
HOURS=$(seq 0 23)
MINUTES=$(seq 0 5 55)

# loop through all years
for YEAR in ${YEARS}
do
	year=${YEAR}
	determine_month_and_day ${YEAR} ${SLURM_ARRAY_TASK_ID}
	month=${array_out[0]}
	day=${array_out[1]}

	# loop through all hours and minutes of this day in $year and unzip data
	if [[ $month != "NULL" ]] && [[ $day != "NULL" ]] # not day 366 of a year with only 365 days
	then
		for HOUR in ${HOURS} # loop through all hours
		do
			if [ ${HOUR} -lt 10 ]
			then
				hour=0${HOUR}
			else
				hour=${HOUR}
			fi
			for MINUTE in ${MINUTES} # loop through all minutes
			do
				if [ ${MINUTE} -lt 10 ]
				then
					minute=0${MINUTE}
				else
					minute=${MINUTE}
				fi
				DATETIME=${year}${month}${day}-${hour}${minute}
				# check if a .gz file is present; if so unzip
				FILE=*"${DATETIME}"*".gz"
				if compgen -G "$FILE" > /dev/null; then
					gunzip $FILE
					echo "Unzipped data for datetime: ${DATETIME}"
				fi
			done
		done
	fi
done