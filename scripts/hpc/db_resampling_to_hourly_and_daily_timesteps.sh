#!/bin/bash
#SBATCH -D /project/quinnlab/dcl3nd/norfolk		 # working directory
#SBATCH -o /project/quinnlab/dcl3nd/norfolk/scripts/script_out_db/job.d1.%j.%N.out	# Name of the output file (eg. myMPI.oJobID)
#SBATCH --ntasks=1				# Number of tasks per serial job (must be 1)
#SBATCH -p standard				# Queue name "standard" (serial)
#SBATCH -A quinnlab_paid				# allocation name
#SBATCH -t 36:00:00				# Run time per serial job (hh:mm:ss)
#SBATCH --array=1-366			# Array of jobs to loop through (366 days)
#SBATCH --mem-per-cpu=30000

# sacct -j <jobid> -o jobid,jobname%20,user,partition,state,start,end,NodeList%60

# d4a_resampling_to_hourly_and_daily_timesteps.sh

module purge
module load anaconda
source activate mrms_processing

# month is 01 and day is ${SLURM_ARRAY_TASK_ID} for first 31 days of the year
if [ ${SLURM_ARRAY_TASK_ID} -le 31 ]
then
	month=01
	if [ ${SLURM_ARRAY_TASK_ID} -lt 10 ]
	then
		day=0${SLURM_ARRAY_TASK_ID}
	else
		day=${SLURM_ARRAY_TASK_ID}
	fi
fi

# all years, hours and minutes to loop through for each day of the year
YEARS=$(seq 2001 2022)

# loop through all years
for YEAR in ${YEARS}
do
	year=${YEAR}
	# figure out month and day of month from day of year
	if [[ ${SLURM_ARRAY_TASK_ID} -gt 31 ]] # day is past January
	then
		if [[ $(( ${YEAR} % 4 )) -gt 0 ]] # not leap year
		then
			if [[ ${SLURM_ARRAY_TASK_ID} -le 31+28 ]] # day is in February
			then
				month=02
				if [[ ${SLURM_ARRAY_TASK_ID}-31 -lt 10 ]]
				then
					day=0$((${SLURM_ARRAY_TASK_ID}-31))
				else
					day=$((${SLURM_ARRAY_TASK_ID}-31))
				fi
			elif [[ ${SLURM_ARRAY_TASK_ID} -le 31+28+31 ]] # day is in March
			then
				month=03
				if [[ ${SLURM_ARRAY_TASK_ID}-31-28 -lt 10 ]]
				then
					day=0$((${SLURM_ARRAY_TASK_ID}-31-28))
				else
					day=$((${SLURM_ARRAY_TASK_ID}-31-28))
				fi
			elif [[ ${SLURM_ARRAY_TASK_ID} -le 31+28+31+30 ]] # day is in April
			then
				month=04
				if [[ ${SLURM_ARRAY_TASK_ID}-31-28-31 -lt 10 ]]
				then
					day=0$((${SLURM_ARRAY_TASK_ID}-31-28-31))
				else
					day=$((${SLURM_ARRAY_TASK_ID}-31-28-31))
				fi
			elif [[ ${SLURM_ARRAY_TASK_ID} -le 31+28+31+30+31 ]] # day is in May
			then
				month=05
				if [[ ${SLURM_ARRAY_TASK_ID}-31-28-31-30 -lt 10 ]]
				then
					day=0$((${SLURM_ARRAY_TASK_ID}-31-28-31-30))
				else
					day=$((${SLURM_ARRAY_TASK_ID}-31-28-31-30))
				fi
			elif [[ ${SLURM_ARRAY_TASK_ID} -le 31+28+31+30+31+30 ]] # day is in June
			then
				month=06
				if [[ ${SLURM_ARRAY_TASK_ID}-31-28-31-30-31 -lt 10 ]]
				then
					day=0$((${SLURM_ARRAY_TASK_ID}-31-28-31-30-31))
				else
					day=$((${SLURM_ARRAY_TASK_ID}-31-28-31-30-31))
				fi
			elif [[ ${SLURM_ARRAY_TASK_ID} -le 31+28+31+30+31+30+31 ]] # day is in July
			then
				month=07
				if [[ ${SLURM_ARRAY_TASK_ID}-31-28-31-30-31-30 -lt 10 ]]
				then
					day=0$((${SLURM_ARRAY_TASK_ID}-31-28-31-30-31-30))
				else
					day=$((${SLURM_ARRAY_TASK_ID}-31-28-31-30-31-30))
				fi
			elif [[ ${SLURM_ARRAY_TASK_ID} -le 31+28+31+30+31+30+31+31 ]] # day is in August
			then
				month=08
				if [[ ${SLURM_ARRAY_TASK_ID}-31-28-31-30-31-30-31 -lt 10 ]]
				then
					day=0$((${SLURM_ARRAY_TASK_ID}-31-28-31-30-31-30-31))
				else
					day=$((${SLURM_ARRAY_TASK_ID}-31-28-31-30-31-30-31))
				fi
			elif [[ ${SLURM_ARRAY_TASK_ID} -le 31+28+31+30+31+30+31+31+30 ]] # day is in September
			then
				month=09
				if [[ ${SLURM_ARRAY_TASK_ID}-31-28-31-30-31-30-31-31 -lt 10 ]]
				then
					day=0$((${SLURM_ARRAY_TASK_ID}-31-28-31-30-31-30-31-31))
				else
					day=$((${SLURM_ARRAY_TASK_ID}-31-28-31-30-31-30-31-31))
				fi
			elif [[ ${SLURM_ARRAY_TASK_ID} -le 31+28+31+30+31+30+31+31+30+31 ]] # day is in October
			then
				month=10
				if [[ ${SLURM_ARRAY_TASK_ID}-31-28-31-30-31-30-31-31-30 -lt 10 ]]
				then
					day=0$((${SLURM_ARRAY_TASK_ID}-31-28-31-30-31-30-31-31-30))
				else
					day=$((${SLURM_ARRAY_TASK_ID}-31-28-31-30-31-30-31-31-30))
				fi
			elif [[ ${SLURM_ARRAY_TASK_ID} -le 31+28+31+30+31+30+31+31+30+31+30 ]] # day is in November
			then
				month=11
				if [[ ${SLURM_ARRAY_TASK_ID}-31-28-31-30-31-30-31-31-30-31 -lt 10 ]]
				then
					day=0$((${SLURM_ARRAY_TASK_ID}-31-28-31-30-31-30-31-31-30-31))
				else
					day=$((${SLURM_ARRAY_TASK_ID}-31-28-31-30-31-30-31-31-30-31))
				fi
			elif [[ ${SLURM_ARRAY_TASK_ID} -le 31+28+31+30+31+30+31+31+30+31+30+31 ]] # day is in December
			then
				month=12
				if [[ ${SLURM_ARRAY_TASK_ID}-31-28-31-30-31-30-31-31-30-31-30 -lt 10 ]]
				then
					day=0$((${SLURM_ARRAY_TASK_ID}-31-28-31-30-31-30-31-31-30-31-30))
				else
					day=$((${SLURM_ARRAY_TASK_ID}-31-28-31-30-31-30-31-31-30-31-30))
				fi
			else # day is 366 but there are only 365 days in the year
				month=NULL
				day=NULL
			fi
		else # leap year
			if [[ ${SLURM_ARRAY_TASK_ID} -le 31+29 ]] # day is in February
			then
				month=02
				if [[ ${SLURM_ARRAY_TASK_ID}-31 -lt 10 ]]
				then
					day=0$((${SLURM_ARRAY_TASK_ID}-31))
				else
					day=$((${SLURM_ARRAY_TASK_ID}-31))
				fi
			elif [[ ${SLURM_ARRAY_TASK_ID} -le 31+29+31 ]] # day is in March
			then
				month=03
				if [[ ${SLURM_ARRAY_TASK_ID}-31-29 -lt 10 ]]
				then
					day=0$((${SLURM_ARRAY_TASK_ID}-31-29))
				else
					day=$((${SLURM_ARRAY_TASK_ID}-31-29))
				fi
			elif [[ ${SLURM_ARRAY_TASK_ID} -le 31+29+31+30 ]] # day is in April
			then
				month=04
				if [[ ${SLURM_ARRAY_TASK_ID}-31-29-31 -lt 10 ]]
				then
					day=0$((${SLURM_ARRAY_TASK_ID}-31-29-31))
				else
					day=$((${SLURM_ARRAY_TASK_ID}-31-29-31))
				fi
			elif [[ ${SLURM_ARRAY_TASK_ID} -le 31+29+31+30+31 ]] # day is in May
			then
				month=05
				if [[ ${SLURM_ARRAY_TASK_ID}-31-29-31-30 -lt 10 ]]
				then
					day=0$((${SLURM_ARRAY_TASK_ID}-31-29-31-30))
				else
					day=$((${SLURM_ARRAY_TASK_ID}-31-29-31-30))
				fi
			elif [[ ${SLURM_ARRAY_TASK_ID} -le 31+29+31+30+31+30 ]] # day is in June
			then
				month=06
				if [[ ${SLURM_ARRAY_TASK_ID}-31-29-31-30-31 -lt 10 ]]
				then
					day=0$((${SLURM_ARRAY_TASK_ID}-31-29-31-30-31))
				else
					day=$((${SLURM_ARRAY_TASK_ID}-31-29-31-30-31))
				fi
			elif [[ ${SLURM_ARRAY_TASK_ID} -le 31+29+31+30+31+30+31 ]] # day is in July
			then
				month=07
				if [[ ${SLURM_ARRAY_TASK_ID}-31-29-31-30-31-30 -lt 10 ]]
				then
					day=0$((${SLURM_ARRAY_TASK_ID}-31-29-31-30-31-30))
				else
					day=$((${SLURM_ARRAY_TASK_ID}-31-29-31-30-31-30))
				fi
			elif [[ ${SLURM_ARRAY_TASK_ID} -le 31+29+31+30+31+30+31+31 ]] # day is in August
			then
				month=08
				if [[ ${SLURM_ARRAY_TASK_ID}-31-29-31-30-31-30-31 -lt 10 ]]
				then
					day=0$((${SLURM_ARRAY_TASK_ID}-31-29-31-30-31-30-31))
				else
					day=$((${SLURM_ARRAY_TASK_ID}-31-29-31-30-31-30-31))
				fi
			elif [[ ${SLURM_ARRAY_TASK_ID} -le 31+29+31+30+31+30+31+31+30 ]] # day is in September
			then
				month=09
				if [[ ${SLURM_ARRAY_TASK_ID}-31-29-31-30-31-30-31-31 -lt 10 ]]
				then
					day=0$((${SLURM_ARRAY_TASK_ID}-31-29-31-30-31-30-31-31))
				else
					day=$((${SLURM_ARRAY_TASK_ID}-31-29-31-30-31-30-31-31))
				fi
			elif [[ ${SLURM_ARRAY_TASK_ID} -le 31+29+31+30+31+30+31+31+30+31 ]] # day is in October
			then
				month=10
				if [[ ${SLURM_ARRAY_TASK_ID}-31-29-31-30-31-30-31-31-30 -lt 10 ]]
				then
					day=0$((${SLURM_ARRAY_TASK_ID}-31-29-31-30-31-30-31-31-30))
				else
					day=$((${SLURM_ARRAY_TASK_ID}-31-29-31-30-31-30-31-31-30))
				fi
			elif [[ ${SLURM_ARRAY_TASK_ID} -le 31+29+31+30+31+30+31+31+30+31+30 ]] # day is in November
			then
				month=11
				if [[ ${SLURM_ARRAY_TASK_ID}-31-29-31-30-31-30-31-31-30-31 -lt 10 ]]
				then
					day=0$((${SLURM_ARRAY_TASK_ID}-31-29-31-30-31-30-31-31-30-31))
				else
					day=$((${SLURM_ARRAY_TASK_ID}-31-29-31-30-31-30-31-31-30-31))
				fi
			elif [[ ${SLURM_ARRAY_TASK_ID} -le 31+29+31+30+31+30+31+31+30+31+30+31 ]] # day is in December
			then
				month=12
				if [[ ${SLURM_ARRAY_TASK_ID}-31-29-31-30-31-30-31-31-30-31-30 -lt 10 ]]
				then
					day=0$((${SLURM_ARRAY_TASK_ID}-31-29-31-30-31-30-31-31-30-31-30))
				else
					day=$((${SLURM_ARRAY_TASK_ID}-31-29-31-30-31-30-31-31-30-31-30))
				fi
			fi
		fi
	fi

 echo "Node ID: $HOSTNAME"
 echo "Slurm Array Task ID: ${SLURM_ARRAY_TASK_ID}"
 echo "Attempting to create hourly and daily netcdfs for ${year}${month}${day}"
 python scripts/_d4a_resampling_to_hourly_and_daily_timesteps.py ${year}${month}${day} "/project/quinnlab/dcl3nd/norfolk/data/mrms_for_rainyday/" "/project/quinnlab/dcl3nd/norfolk/data/_scratch_zarrs/" "/project/quinnlab/dcl3nd/norfolk/data/mrms_for_rainyday_hourly/" "/project/quinnlab/dcl3nd/norfolk/data/mrms_for_rainyday_daily/"
done

