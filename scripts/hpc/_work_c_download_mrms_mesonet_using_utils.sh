#!/bin/bash
#SBATCH -D /project/quinnlab/dcl3nd/norfolk/data/raw_data/mrms_grib_mesonet	 # working directory
#SBATCH -o /project/quinnlab/dcl3nd/norfolk/scripts/script_out_c/job.%j.%N.out	# Name of the output file (eg. myMPI.oJobID)
#SBATCH --ntasks=1				# Number of tasks per serial job (must be 1)
#SBATCH -p standard				# Queue name "standard" (serial)
#SBATCH -A quinnlab_paid				# allocation name
#SBATCH -t 48:00:00				# Run time per serial job (hh:mm:ss)
#SBATCH --array=1-366%20			# Array of jobs to loop through (366 days)

source _work__utils.sh # this contains the code to return the month and day given a year and array id

# all years, hours and minutes to loop through for each day of the year
YEARS=$(seq 2015 2022)
HOURS=$(seq 0 23)
MINUTES=$(seq 0 2 58)

# loop through all years
for YEAR in ${YEARS}
do
	year=${YEAR}
	determine_month_and_day ${YEAR} ${SLURM_ARRAY_TASK_ID}
	month=${array_out[0]}
	day=${array_out[1]}
	
	# loop through all hours and minutes of this day in $year and download and unzip data
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
				# download and unzip .grib.gz if it's not already there
				## check if the unzipped file already exits
				FILE="PrecipRate_00.00_${year}${month}${day}-${hour}${minute}00.grib2"
				echo "File being processed: $FILE"
				# Download and unzip the file only if the .grib file doesn't exist
				# source: https://stackoverflow.com/questions/6363441/check-if-a-file-exists-with-a-wildcard-in-a-shell-script
				if compgen -G "$FILE" > /dev/null; then
					echo "File already exists, skipping the downloading and unzipping..."
				else
					echo "File does not exist! Downloading and unzipping..."
					wget -q https://mtarchive.geol.iastate.edu/${year}/${month}/${day}/mrms/ncep/PrecipRate/PrecipRate_00.00_${year}${month}${day}-${hour}${minute}00.grib2.gz
					gunzip PrecipRate_00.00_${year}${month}${day}-${hour}${minute}00.grib2.gz
					rm *${year}${month}${day}-${hour}${minute}*.gz*
					echo "Downloaded and unzipped file."
				fi
				# check if a .gz file is present; if so remove it
				FILE=*${year}${month}${day}-${hour}${minute}*.gz*
				if compgen -G "$FILE" > /dev/null; then
					echo "A .gz file is present; removing..."
					rm *${year}${month}${day}-${hour}${minute}*.gz*
				else
					echo "No .gz file present. Script is complete."
				fi
				echo "done!"
			done
		done
	fi
done