#!/bin/bash
# Function to determine the first job in the array
is_first_array_job() {
    MASTER_JOB_ID=${SLURM_ARRAY_JOB_ID:-$SLURM_JOB_ID}
    ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}

    # Query the job array IDs using SLURM environment variables (case-insensitive parsing)
    local first_job=$(scontrol show job ${MASTER_JOB_ID} | grep -iPo 'ArrayTaskId=\K[0-9]+' | sort -n | head -1)

    # Handle case where the first job ID is not found
    if [[ -z "$first_job" ]]; then
        echo "Warning: Unable to retrieve the first job ID. Defaulting to task ID 1."
        first_job=1
    fi

    echo "First job in the array: $first_job"

    # Check if the current task ID matches the first job in the array
    if [ "$ARRAY_TASK_ID" -eq "$first_job" ]; then
        return 0  # True: this is the first job in the array
    else
        return 1  # False: this is not the first job
    fi
}

# Only the first job in the array performs the archiving
archive_previous_script_outfiles() {
    if is_first_array_job; then
        mkdir -p ${dir_outs}${SLURM_JOB_NAME}/_archive
        echo "Archiving outputs... Only task ${ARRAY_TASK_ID} is responsible for this."

        find ${dir_outs}${SLURM_JOB_NAME} -maxdepth 1 -name "*.out" ! -name "${MASTER_JOB_ID}*.out" \
            -exec mv {} ${dir_outs}${SLURM_JOB_NAME}/_archive/ \; 2>/dev/null || true
    else
        echo "Skipping archiving for task ${ARRAY_TASK_ID}. Not the first job in the array."
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
