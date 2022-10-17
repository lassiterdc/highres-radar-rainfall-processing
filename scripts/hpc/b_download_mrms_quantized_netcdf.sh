#!/bin/bash
#SBATCH -D /project/quinnlab/dcl3nd/norfolk/data/raw_data/mrms_nc_quant	 # working directory
#SBATCH -o /project/quinnlab/dcl3nd/norfolk/scripts/script_out_b/job.%j.%N.out	# Name of the output file (eg. myMPI.oJobID)
#SBATCH --ntasks=1				# Number of tasks per serial job (must be 1)
#SBATCH -p standard				# Queue name "standard" (serial)
#SBATCH -A quinnlab_paid				# allocation name
#SBATCH -t 48:00:00				# Run time per serial job (hh:mm:ss)
#SBATCH --array=1-366%20			# Array of jobs to loop through (366 days)

# Once script has run, delete the empty files (this looped through every minute since portions are at 2 minutes and portions are at 5 minutes)
# https://www.cyberciti.biz/faq/howto-find-delete-empty-directories-files-in-unix-linux/
# find /path/to/dir -empty -type f -delete

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

# all years
YEARS=$(seq 2012 2014)

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

  # loop through all hours and minutes of this day in $year and download and unzip data
  HOURS=$(seq 0 23)
  MINUTES=$(seq 0 1 59)
  url1='https://mesonet.agron.iastate.edu/cgi-bin/request/raster2netcdf.py?dstr='
  url2='&prod=mrms_a2m'
  # source for URL: https://mesonet.agron.iastate.edu/GIS/rasters.php?rid=3
	if [[ $month != "NULL" ]] && [[ $day != "NULL" ]] # not day 366 of a year with only 365 days
	then   
    # loop through all hours and minutes
		for HOUR in ${HOURS}
		do
			if [ ${HOUR} -lt 10 ]
			then
				hour=0${HOUR}
			else
				hour=${HOUR}
			fi
			for MINUTE in ${MINUTES} # minutes
			do
				if [ ${MINUTE} -lt 10 ]
				then
					minute=0${MINUTE}
				else
					minute=${MINUTE}
				fi
				# download file only if it is not there
				## check if the file already exits
				datetime=${year}${month}${day}${hour}${minute}
				FILE="a2m_${datetime}.nc"
				echo "File being processed: $FILE"
				# Download the file only if it doesn't exist
				# source: https://stackoverflow.com/questions/6363441/check-if-a-file-exists-with-a-wildcard-in-a-shell-script
				if compgen -G "$FILE" > /dev/null; then
					echo "File already exists, skipping the downloading..."
				else
					echo "File does not exist! Downloading..."
					full_url=${url1}${datetime}${url2}
					wget -q -O a2m_${datetime}.nc ${full_url}
					echo "downloaded data for ${year}${month}${day}-${hour}${minute}"
				fi
				echo "done!"
				find . -name "$FILE" -type f -size -2k -delete # delete file if it is below 2kilobytes
			done
		done
	fi
done