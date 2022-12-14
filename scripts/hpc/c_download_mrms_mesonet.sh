#!/bin/bash
#SBATCH -o _script_outputs/%x/%A_%a_%N.out
#SBATCH -e _script_errors/%x/%A_%a_%N.out
#SBATCH --ntasks=1				# Number of tasks per serial job (must be 1)
#SBATCH -p standard				# Queue name "standard" (serial)
#SBATCH -A quinnlab_paid				# allocation name
#SBATCH -t 72:00:00				# Run time per serial job (hh:mm:ss)
#SBATCH --array=295,296,297,298,299,300,301,302,303,304,305,306,307,308,309,312,313,314,315,316,317,318,319,320
#SBATCH --mail-user=dcl3nd@virginia.edu          # address for email notification
#SBATCH --mail-type=ALL   
#SBATCH --exclude=udc-aw29-25b,udc-an33-5c0,udc-an33-7c1,udc-aw29-19b,udc-an33-11c1,udc-aw34-3c0,udc-ba26-34c1,udc-aw34-4c0,udc-ba25-32c1,udc-aw29-23a,udc-aw34-19c0,udc-aw34-11c1,udc-aw34-3c1						
# add %20 after #SBATCH --array=1-366 if this is the first passthrough

source __utils.sh
source __directories.sh
#confirm working directory exists
mkdir -p ${assar_dirs[repo]}${assar_dirs[raw_mrms]}
# move to working directory
cd ${assar_dirs[repo]}${assar_dirs[raw_mrms]}

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
				# start timer at 0
				SECONDS=0
				DATETIME=${year}${month}${day}-${hour}${minute}
				FILE=*"${DATETIME}"*".grib2"
				# if the .grib file does not exist, download it
				if ! compgen -G "$FILE" > /dev/null; then
					wget -q -c https://mtarchive.geol.iastate.edu/${year}/${month}/${day}/mrms/ncep/PrecipRate/PrecipRate_00.00_${DATETIME}00.grib2.gz
					downloaded="was"
					# echo "Downloaded data for datetime: ${DATETIME}"
				else
					# echo ".grib file already exists for datetime: ${DATETIME}"
					downloaded="was not"
				fi
				# if a .gz file exist, unzip it
				FILE=*"${DATETIME}"*".gz"
				if compgen -G "$FILE" > /dev/null; then
					gunzip $FILE
					# echo "Unzipped .gz file for datetime: ${DATETIME}"
					unzipped="was"
				else
					# echo "No .gz file present for datetime: ${DATETIME}"
					unzipped="was not"
				fi
				duration=$SECONDS
				echo "Processed datetime ${DATETIME}; Time elapsed: $(($duration / 60)) minutes and $(($duration % 60)) seconds; file $downloaded downloaded and $unzipped unzipped."
			done
		done
	fi
done