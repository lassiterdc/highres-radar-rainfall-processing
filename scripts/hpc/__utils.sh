#!/bin/bash

move_out_and_error_files () {
	# OUTFILE=$1
	OUTFILE="${0##*/}.out"
	OUTFOLDER=${assar_dirs[repo]}${assar_dirs[hpc_outputs]}
	ERRORFOLDER=${assar_dirs[repo]}${assar_dirs[hpc_errors]}

	ARCHIVEFOLDER="_archive/"

	if compgen -G "$OUTFOLDER$OUTFILE" > /dev/null; then
	# echo "File already exists, moving to archive..."
	ARCHIVEFOLDER_FULLPATH = $OUTFOLDER$ARCHIVEFOLDER
	mkdir -p $ARCHIVEFOLDER_FULLPATH
	mv $OUTFOLDER$OUTFILE $ARCHIVEFOLDER_FULLPATH$OUTFILE
	fi

	if compgen -G "$ERRORFOLDER$OUTFILE" > /dev/null; then
	# echo "File already exists, moving to archive..."
	ARCHIVEFOLDER_FULLPATH = $ERRORFOLDER$ARCHIVEFOLDER
	mkdir -p $ARCHIVEFOLDER_FULLPATH
	mv $ERRORFOLDER$OUTFILE $ARCHIVEFOLDER_FULLPATH$OUTFILE
	fi
}

# The following function takes a four-digit year and SLURM_ARRAY_TASK_ID and returns
# an array named array_out whose 0th index is the month 1st index is the day taking leap years into account
determine_month_and_day () {
	YEAR=$1 # first argument
	TASK_ID=$2 # second argument
	# MONTH is 01 and DAY is ${SLURM_ARRAY_TASK_ID} for first 31 DAYs of the YEAR
	if [ ${TASK_ID} -le 31 ]
	then
		MONTH=01
		if [ ${TASK_ID} -lt 10 ]
		then
			DAY=0${TASK_ID}
		else
			DAY=${TASK_ID}
		fi
	fi

	# figure out MONTH and DAY of MONTH from DAY of YEAR
	if [[ ${TASK_ID} -gt 31 ]] # DAY is past January
	then
		if [[ $(( ${YEAR} % 4 )) -gt 0 ]] # not leap YEAR
		then
			if [[ ${TASK_ID} -le 31+28 ]] # DAY is in February
			then
				MONTH=02
				if [[ ${TASK_ID}-31 -lt 10 ]]
				then
					DAY=0$((${TASK_ID}-31))
				else
					DAY=$((${TASK_ID}-31))
				fi
			elif [[ ${TASK_ID} -le 31+28+31 ]] # DAY is in March
			then
				MONTH=03
				if [[ ${TASK_ID}-31-28 -lt 10 ]]
				then
					DAY=0$((${TASK_ID}-31-28))
				else
					DAY=$((${TASK_ID}-31-28))
				fi
			elif [[ ${TASK_ID} -le 31+28+31+30 ]] # DAY is in April
			then
				MONTH=04
				if [[ ${TASK_ID}-31-28-31 -lt 10 ]]
				then
					DAY=0$((${TASK_ID}-31-28-31))
				else
					DAY=$((${TASK_ID}-31-28-31))
				fi
			elif [[ ${TASK_ID} -le 31+28+31+30+31 ]] # DAY is in May
			then
				MONTH=05
				if [[ ${TASK_ID}-31-28-31-30 -lt 10 ]]
				then
					DAY=0$((${TASK_ID}-31-28-31-30))
				else
					DAY=$((${TASK_ID}-31-28-31-30))
				fi
			elif [[ ${TASK_ID} -le 31+28+31+30+31+30 ]] # DAY is in June
			then
				MONTH=06
				if [[ ${TASK_ID}-31-28-31-30-31 -lt 10 ]]
				then
					DAY=0$((${TASK_ID}-31-28-31-30-31))
				else
					DAY=$((${TASK_ID}-31-28-31-30-31))
				fi
			elif [[ ${TASK_ID} -le 31+28+31+30+31+30+31 ]] # DAY is in July
			then
				MONTH=07
				if [[ ${TASK_ID}-31-28-31-30-31-30 -lt 10 ]]
				then
					DAY=0$((${TASK_ID}-31-28-31-30-31-30))
				else
					DAY=$((${TASK_ID}-31-28-31-30-31-30))
				fi
			elif [[ ${TASK_ID} -le 31+28+31+30+31+30+31+31 ]] # DAY is in August
			then
				MONTH=08
				if [[ ${TASK_ID}-31-28-31-30-31-30-31 -lt 10 ]]
				then
					DAY=0$((${TASK_ID}-31-28-31-30-31-30-31))
				else
					DAY=$((${TASK_ID}-31-28-31-30-31-30-31))
				fi
			elif [[ ${TASK_ID} -le 31+28+31+30+31+30+31+31+30 ]] # DAY is in September
			then
				MONTH=09
				if [[ ${TASK_ID}-31-28-31-30-31-30-31-31 -lt 10 ]]
				then
					DAY=0$((${TASK_ID}-31-28-31-30-31-30-31-31))
				else
					DAY=$((${TASK_ID}-31-28-31-30-31-30-31-31))
				fi
			elif [[ ${TASK_ID} -le 31+28+31+30+31+30+31+31+30+31 ]] # DAY is in October
			then
				MONTH=10
				if [[ ${TASK_ID}-31-28-31-30-31-30-31-31-30 -lt 10 ]]
				then
					DAY=0$((${TASK_ID}-31-28-31-30-31-30-31-31-30))
				else
					DAY=$((${TASK_ID}-31-28-31-30-31-30-31-31-30))
				fi
			elif [[ ${TASK_ID} -le 31+28+31+30+31+30+31+31+30+31+30 ]] # DAY is in November
			then
				MONTH=11
				if [[ ${TASK_ID}-31-28-31-30-31-30-31-31-30-31 -lt 10 ]]
				then
					DAY=0$((${TASK_ID}-31-28-31-30-31-30-31-31-30-31))
				else
					DAY=$((${TASK_ID}-31-28-31-30-31-30-31-31-30-31))
				fi
			elif [[ ${TASK_ID} -le 31+28+31+30+31+30+31+31+30+31+30+31 ]] # DAY is in December
			then
				MONTH=12
				if [[ ${TASK_ID}-31-28-31-30-31-30-31-31-30-31-30 -lt 10 ]]
				then
					DAY=0$((${TASK_ID}-31-28-31-30-31-30-31-31-30-31-30))
				else
					DAY=$((${TASK_ID}-31-28-31-30-31-30-31-31-30-31-30))
				fi
			else # DAY is 366 but there are only 365 DAYs in the YEAR
				MONTH=NULL
				DAY=NULL
			fi
		else # leap YEAR
			if [[ ${TASK_ID} -le 31+29 ]] # DAY is in February
			then
				MONTH=02
				if [[ ${TASK_ID}-31 -lt 10 ]]
				then
					DAY=0$((${TASK_ID}-31))
				else
					DAY=$((${TASK_ID}-31))
				fi
			elif [[ ${TASK_ID} -le 31+29+31 ]] # DAY is in March
			then
				MONTH=03
				if [[ ${TASK_ID}-31-29 -lt 10 ]]
				then
					DAY=0$((${TASK_ID}-31-29))
				else
					DAY=$((${TASK_ID}-31-29))
				fi
			elif [[ ${TASK_ID} -le 31+29+31+30 ]] # DAY is in April
			then
				MONTH=04
				if [[ ${TASK_ID}-31-29-31 -lt 10 ]]
				then
					DAY=0$((${TASK_ID}-31-29-31))
				else
					DAY=$((${TASK_ID}-31-29-31))
				fi
			elif [[ ${TASK_ID} -le 31+29+31+30+31 ]] # DAY is in May
			then
				MONTH=05
				if [[ ${TASK_ID}-31-29-31-30 -lt 10 ]]
				then
					DAY=0$((${TASK_ID}-31-29-31-30))
				else
					DAY=$((${TASK_ID}-31-29-31-30))
				fi
			elif [[ ${TASK_ID} -le 31+29+31+30+31+30 ]] # DAY is in June
			then
				MONTH=06
				if [[ ${TASK_ID}-31-29-31-30-31 -lt 10 ]]
				then
					DAY=0$((${TASK_ID}-31-29-31-30-31))
				else
					DAY=$((${TASK_ID}-31-29-31-30-31))
				fi
			elif [[ ${TASK_ID} -le 31+29+31+30+31+30+31 ]] # DAY is in July
			then
				MONTH=07
				if [[ ${TASK_ID}-31-29-31-30-31-30 -lt 10 ]]
				then
					DAY=0$((${TASK_ID}-31-29-31-30-31-30))
				else
					DAY=$((${TASK_ID}-31-29-31-30-31-30))
				fi
			elif [[ ${TASK_ID} -le 31+29+31+30+31+30+31+31 ]] # DAY is in August
			then
				MONTH=08
				if [[ ${TASK_ID}-31-29-31-30-31-30-31 -lt 10 ]]
				then
					DAY=0$((${TASK_ID}-31-29-31-30-31-30-31))
				else
					DAY=$((${TASK_ID}-31-29-31-30-31-30-31))
				fi
			elif [[ ${TASK_ID} -le 31+29+31+30+31+30+31+31+30 ]] # DAY is in September
			then
				MONTH=09
				if [[ ${TASK_ID}-31-29-31-30-31-30-31-31 -lt 10 ]]
				then
					DAY=0$((${TASK_ID}-31-29-31-30-31-30-31-31))
				else
					DAY=$((${TASK_ID}-31-29-31-30-31-30-31-31))
				fi
			elif [[ ${TASK_ID} -le 31+29+31+30+31+30+31+31+30+31 ]] # DAY is in October
			then
				MONTH=10
				if [[ ${TASK_ID}-31-29-31-30-31-30-31-31-30 -lt 10 ]]
				then
					DAY=0$((${TASK_ID}-31-29-31-30-31-30-31-31-30))
				else
					DAY=$((${TASK_ID}-31-29-31-30-31-30-31-31-30))
				fi
			elif [[ ${TASK_ID} -le 31+29+31+30+31+30+31+31+30+31+30 ]] # DAY is in November
			then
				MONTH=11
				if [[ ${TASK_ID}-31-29-31-30-31-30-31-31-30-31 -lt 10 ]]
				then
					DAY=0$((${TASK_ID}-31-29-31-30-31-30-31-31-30-31))
				else
					DAY=$((${TASK_ID}-31-29-31-30-31-30-31-31-30-31))
				fi
			elif [[ ${TASK_ID} -le 31+29+31+30+31+30+31+31+30+31+30+31 ]] # DAY is in December
			then
				MONTH=12
				if [[ ${TASK_ID}-31-29-31-30-31-30-31-31-30-31-30 -lt 10 ]]
				then
					DAY=0$((${TASK_ID}-31-29-31-30-31-30-31-31-30-31-30))
				else
					DAY=$((${TASK_ID}-31-29-31-30-31-30-31-31-30-31-30))
				fi
			fi
		fi
	fi
	# initiate the output array
	array_out[0]=$MONTH
	array_out[1]=$DAY
	# echo "Month: $MONTH"
	# echo "Day: $DAY"
	# echo ${array_out[@]}
	# echo ${array_out[@]}
}
