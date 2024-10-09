#!/bin/bash
#SBATCH -o _script_outputs/%x/%A_%a_%N.out
#SBATCH -e _script_errors/%x/%A_%a_%N.out
#SBATCH --ntasks=1				# Number of tasks per serial job (must be 1)
#SBATCH -p standard				# Queue name "standard" (serial)
#SBATCH -A quinnlab				# allocation name
#SBATCH -t 168:00:00				# Run time per serial job (hh:mm:ss)
#SBATCH --array=1-366%20 # 295,296,297,298,299,300,301,302,303,304,305,306,307,308,309,312,313,314,315,316,317,318,319,320
#SBATCH --mail-user=dcl3nd@virginia.edu          # address for email notification
#SBATCH --mail-type=ALL   
# SBATCH --exclude=udc-aw29-25b,udc-an33-5c0,udc-an33-7c1,udc-aw29-19b,udc-an33-11c1,udc-aw34-3c0,udc-ba26-34c1,udc-aw34-4c0,udc-ba25-32c1,udc-aw29-23a,udc-aw34-19c0,udc-aw34-11c1,udc-aw34-3c1						
# add %20 after #SBATCH --array=1-366 if this is the first passthrough

# ijob -c 1 -A quinnlab -p standard --time=0-08:00:00

# cd /project/quinnlab/dcl3nd/norfolk/highres-radar-rainfall-processing/scripts/hpc
# rm /project/quinnlab/dcl3nd/norfolk/highres-radar-rainfall-processing/scripts/hpc/c_download_mrms_mesonet.sh
# git pull

source __utils.sh
source __directories.sh
#confirm working directory exists
mkdir -p ${assar_dirs[raw_mrms]}
# move to working directory
cd ${assar_dirs[raw_mrms]}

# all years, hours and minutes to loop through for each day of the year
YEARS=$(seq 2015 2024)

# loop through all years
for YEAR in ${YEARS}; do
year=${YEAR}
determine_month_and_day ${YEAR} ${SLURM_ARRAY_TASK_ID}
month=${array_out[0]}
day=${array_out[1]}
# start timer at 0
SECONDS=0
# return a list of all files for this day and write to a text file
url="https://mtarchive.geol.iastate.edu/${year}/${month}/${day}/mrms/ncep/PrecipRate/"
echo "$url"
wget --quiet -r --no-parent -l1 -nd -A "*.grib2.gz" "$url" -O directory_listing_${year}${month}${day}.html
grep -a -oP 'href="\K[^"]+\.grib2\.gz' directory_listing_${year}${month}${day}.html > file_list_${year}${month}${day}.txt
# sort contents and make sure all lines are unique
sort file_list_${year}${month}${day}.txt | uniq > sorted_file_list_${year}${month}${day}.txt
# Loop through the file list; if they exist locally, skip download; otherwise, download
if [[ ! -f "sorted_file_list_${year}${month}${day}.txt" ]]; then
    echo "Error: File sorted_file_list_${year}${month}${day}.txt not found!" >&2
    continue
fi
counter=0 # for troubleshooting
while IFS= read -r line; do
    counter=$((counter + 1))
    # if [ "$counter" -ge 3 ]; then
    #     break
    # fi
    # Clean file name by removing leading/trailing spaces or newlines
    file=$(echo "$line" | xargs | tr -d '\n')

    # Filter out lines that contain "Binary file" or non-regular lines
    if [[ "$file" == *"Binary file"* ]]; then
        # echo "Skipping binary file line: ${file}"
        continue
    fi

    # Match .grib2.gz files
    if [[ "$file" == *".grib2.gz"* ]]; then
        load_url="${url}${file}"
        local_filename="${assar_dirs[raw_mrms]}${file%.gz}"  # unzipped filename, remove .gz extension
        # Debugging: Show the constructed URL and filename
        # echo "Constructed download URL: '$load_url'"
        # echo "Filename: '$local_filename'"
        # Check if file exists, quoting the variables
        if [[ ! -f "${local_filename}" ]]; then
            # Print the command that will be executed
            # echo "File does not exist locally. Downloading and unzipping..."
            # echo "Running command: wget -q -O - '$load_url' | gunzip -c > '${local_filename}'"

            # Combine download and unzip, ensuring variables are quoted
            if ! wget -q -O - "$load_url" | gunzip -c > "${local_filename}"; then
                echo "Error: Failed to download or unzip '${load_url}'"
                continue  # Skip to the next file if this one fails
            fi
            downloaded="was downloaded and unzipped"
        else
            # echo "File ${local_filename} already exists locally."
            downloaded="was not downloaded because it already exists locally"
        fi
        echo "'${file}' $downloaded."
    fi
done < sorted_file_list_${year}${month}${day}.txt
duration=$SECONDS
echo "Processed date ${month}/${day}/${year}; Time elapsed: $(($duration / 60)) minutes and $(($duration % 60)) seconds"
done