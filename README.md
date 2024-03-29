# highres-radar-rainfall-processing
This repository contains code to download, process, and quality check the [NOAA's Surface Precipitation Rate product](https://vlab.noaa.gov/awips-reference?p_p_id=consolidatedreferencesearch&p_p_lifecycle=0&p_p_state=maximized&p_p_mode=view&_consolidatedreferencesearch_mvcPath=%2Fsearch%2Fview_content.jsp&_consolidatedreferencesearch_redirectBack=https%3A%2F%2Fvlab.noaa.gov%3A443%2Fawips-reference%3Fp_p_id%3Dconsolidatedreferencesearch%26p_p_lifecycle%3D0%26p_p_state%3Dnormal%26p_p_mode%3Dview%26keywords%3Dprecipitation&_consolidatedreferencesearch_assetEntryId=669245&_consolidatedreferencesearch_type=content&_consolidatedreferencesearch_groupId=96675) for use with
[RainyDay](https://github.com/danielbwright/RainyDay) which performs stochastic storm transposition using gridded rainfall products. The repository also contains code for comparing the radar rainfall data with a gage network.

# Data access:
- [Dropbox folder](https://www.dropbox.com/sh/yql8paow2k07aa6/AAA9cIrZzMZnofBTVBBCj-sNa?dl=0) with a sample of MRMS data from 09/22/2009 and 03/06/2019. The temporal resolution of pre-2012 data is 5-minute and post-2012 is 2-minute.
- The full dataset can be made available upon request but sharing requires a Globus subscription

# Key processing steps
- Downloading and unzipping .grib and netcdf files from multiple sources
- Combining the files into daily netcdf files for usability and transferability. Steps include:
    - Quantify the number of missing timesteps each day
    - Ensure consistent variable, dimension, coordinate, and attribute names
    - Verify units and ensure they are clearly documented within the netcdf
    - Other helpful attributes are added the contain information regarding:
        - The source of the data
        - Number of missing timesteps for that day
        - A `warnings` attribute that documents:
            - Changes in the grid coordinates from timestep to timestep
            - Whether the data for that day conform to a standard timestep (2 or 5 minutes)
            - Number of missing timesteps
            - Whether the original grid coordinates of the raw data are identical to the grid coordinates of the most recent raw data file
        - A `problems` attribute that should be blank for most dataset but contains detailed descriptions of potential issues identified in the script like mis-assigned units
    - Negative values used to represent missing or fill values are all assigned a value of `nan`
    - Rainfall intensities of above 9,000 mm/hr are converted to `nan`
    - Fill missing timesteps with `nan` values
    - Ensure that the latitude and longitude coordinates align with the most recent timestep downloaded
        - Where there are differences, assign the most recent latitude and longitude coordinates to the file and log the shift and the original coordinates as an attribute
    - Shift the latitude and longitude coordinates so they represent the upper left corner of each gridcell

# Requirements
- Because the original data has a file for each 2-minute timestep, adding up to millions of files, an HPC was needed to download and process this data
- All in all, it took days of processing time

# How to use:
- If you are interested in using one of the data products, the easiest thing would be to reach out to me and we can figure out the best way to transfer ~3.5TB of data
- If you would like to reproduce the process to download and process the data yourself:
    1. Clone the repository
    2. In `/scripts/hpc/__directories.sh`, change the `assar_dirs[repo]` variable to match the filepath of the cloned repository
    3. In  `/scripts/hpc/`, run the shell scripts in alphabetical order starting at `aa_...` through at least `da_...` to produce the primary output, `mrms_nc_preciprate_fullres_dailyfiles` (see [Outputs](#outputs) for more description)

# Data Inputs
- [The MRMS precipitation rate product ](https://vlab.noaa.gov/awips-reference?p_p_id=consolidatedreferencesearch&p_p_lifecycle=0&p_p_state=maximized&p_p_mode=view&_consolidatedreferencesearch_mvcPath=%2Fsearch%2Fview_content.jsp&_consolidatedreferencesearch_redirectBack=https%3A%2F%2Fvlab.noaa.gov%3A443%2Fawips-reference%3Fp_p_id%3Dconsolidatedreferencesearch%26p_p_lifecycle%3D0%26p_p_state%3Dnormal%26p_p_mode%3Dview%26keywords%3Dprecipitation&_consolidatedreferencesearch_assetEntryId=669245&_consolidatedreferencesearch_type=content&_consolidatedreferencesearch_groupId=96675) was drawn from two sources:
  - 2001-2011:
    - Dr. Jian Zhang And Dr. Jonathan Gourley. (2018). Multi-Radar Multi-Sensor Precipitation Reanalysis (Version 1.0). Open Commons Consortium Environmental Data Commons. https://doi.org/10.25638/EDC.PRECIP.0001
  - 2015-present: [Iowa Environmental Mesonet](https://mesonet.agron.iastate.edu/)
    - Can be accessed through links like this: https://mtarchive.geol.iastate.edu/YYYY/MM/DD/mrms/ncep/PrecipRate/ (e.g., https://mtarchive.geol.iastate.edu/2020/09/09/mrms/ncep/PrecipRate/)

# Outputs
- Datasets
    - The format is `[source]_[type]_[data]_[timestep]_[file aggregation]_[other qualifiers]`
    - CONUS
        - `mrms_nc_preciprate_fullres_dailyfiles`
            - main output; contains a daily netcdf with the full resolution MRMS data
        - `mrms_nc_preciprate_hourly_dailyfiles`
            - daily netcdfs aggregated to an hourly timestep
        - `mrms_nc_preciprate_daily_dailyfiles`
            - daily netcdfs aggregated to a daily timestep
        - `mrms_nc_preciprate_daily_yearlyfiles`
            - yearly netcdfs aggregated to daily timestep
        - `mrms_nc_preciprate_fullres_yearlyfiles_zipped`
            - containis yearly zip folders containing the main output, `mrms_nc_preciprate_fullres_dailyfiles`
        - `mrms_nc_preciprate_yearly_singlefile`
            - single netcdf of yearly rainfall totals
    - Spanning user defined gage network
        - `mrms_nc_preciprate_fullres_yearlyfiles_subset-over-gage-network`
            - yearly netcdfs of full resolution MRMS data subset as a rectangle over over the gage network plus a 5-gridcell buffer on all 4 sides
        - `mrms_csv_preciprate_fullres_yearlyfiles_subset-over-gage-network`
            - yearly csv files of full resolution MRMS at gridcells overlapping gages
- Plots
    - `h_annual_statistics` 
        - Annual rainfall totals and annual rainfall anomolies (difference from mean)
    - Plots of rain gage MRMS comparisons

# Plots
Note: for both plots, `robust=True`, meaning the 2nd and 98th percentiles were used for determining the color bar.

## Annual Rainfall Totals
![annual totals](https://github.com/lassiterdc/highres-radar-rainfall-processing/blob/main/plots/h_annual_statistics/all_years.png?raw=true)

## Annual Rainfall Anomolies (difference from the mean)
![annual anaomolies](https://github.com/lassiterdc/highres-radar-rainfall-processing/blob/main/plots/h_annual_statistics/all_years_anomolies_rel_to_avg_excluding_2012-2014.png?raw=true)

## Gage vs. MRMS Data Comparisons
### Aggregate comparisons
![hist1](https://github.com/lassiterdc/highres-radar-rainfall-processing/blob/main/plots/b_visualizations_of_gage_vs_mrms_data/b_mrms_minus_gage_event_totals_histogram.png?raw=true)

![hist2](https://github.com/lassiterdc/highres-radar-rainfall-processing/blob/main/plots/b_visualizations_of_gage_vs_mrms_data/b_mrms_minus_gage_no_zeros_time_series_histogram.png?raw=true)

![hist3](https://github.com/lassiterdc/highres-radar-rainfall-processing/blob/main/plots/b_visualizations_of_gage_vs_mrms_data/b_mrms_minus_gage_time_series_histogram.png?raw=true)

### Event-based comparisons
#### Time Series Comparisons
![gagevsmrms1](https://github.com/lassiterdc/highres-radar-rainfall-processing/blob/main/plots/b_visualizations_of_gage_vs_mrms_data/MMPS-036/10_e_id_25.png?raw=true)

![gagevsmrms2](https://github.com/lassiterdc/highres-radar-rainfall-processing/blob/main/plots/b_visualizations_of_gage_vs_mrms_data/MMPS-036/14_e_id_151.png?raw=true)

![gagevsmrms3](https://github.com/lassiterdc/highres-radar-rainfall-processing/blob/main/plots/b_visualizations_of_gage_vs_mrms_data/MMPS-036/1_e_id_107.png?raw=true)

#### Spatially Distributed Comparisons
![spatialcomp1](https://github.com/lassiterdc/highres-radar-rainfall-processing/blob/main/plots/b_visualizations_of_gage_vs_mrms_data/mrms/1_e_id_107.png?raw=true)

![spatialcomp2](https://github.com/lassiterdc/highres-radar-rainfall-processing/blob/main/plots/b_visualizations_of_gage_vs_mrms_data/mrms/5_e_id_74.png?raw=true)

![spatialcomp3](https://github.com/lassiterdc/highres-radar-rainfall-processing/blob/main/plots/b_visualizations_of_gage_vs_mrms_data/mrms/10_e_id_25.png?raw=true)

# Information from personal correspondence:
- Jian Zhang from the National Severe Storms Laboratory
    - There have been changes to the algorithms for computing surface precipitation rates from radar reflectivity readings. From 1/1/2014 through 9/30/2020, the precipitation product is based on version 10.5-11.5 (https://doi.org/10.1175/JHM-D-19-0194.1); from 10/1/2020 to the present, the product is based on version version 12-12.1 (https://doi.org/10.1175/JHM-D-19-0194.1)
    - All timesteps are in UTC
    - There should be no IP issues in making my reprocessed MRMS data publically available as long as the original sources are cited
    - The 2001-2011 and 2015-current data have coordinates representing the center of each gridcell
    - The center of the northwestern gridcell for these data should be at 54.995N,129.995W and the center of the southeastern corner grid cell should be at 20.005N, 60.005W
- Daryl Herzmann at the Iowa Environmental Mesonet
  - Data from 2012-August of 2013 are indeed missing
  - The quantized PNGs from August 2013 through 2014 are the only records for that time period.

# Other notes:
- Data seems to be completely missing from the Iowa Environmental Mesonet from 2012 to August of 2013
- From August of 2013 through 2014, a quantized PNG product is available which can be downloaded as a netcdf, but based on annual totals, the process to quantize the data seems to truncate heavy rainfall resulting in underestimation; therefore, data from 2013 and 2014 are omitted from the analysis. However, the code to download and process this data is present in this repository.
  - For information about the quantized PNGs, see this webpage: https://mesonet.agron.iastate.edu/GIS/rasters.php?rid=3
  - download links like this were used to download the quantized PNGs as netcdfs https://mesonet.agron.iastate.edu/cgi-bin/request/raster2netcdf.py?dstr=%Y%m%d%H%M&prod=mrms_a2m (e.g., https://mesonet.agron.iastate.edu/cgi-bin/request/raster2netcdf.py?prod=mrms_a2m&dstr=201710250000)
