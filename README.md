# highres-radar-rainfall-processing
This repository contains code to download, process, and quality check the [NOAA's Surface Precipitation Rate product](https://vlab.noaa.gov/web/wdtd/-/surface-precipitation-rate-spr-?selectedFolder=9234881) for use with
[RainyDay](https://github.com/danielbwright/RainyDay) which performs stochastic storm transposition using gridded rainfall products. The repository also contains code for comparing the radar rainfall data with a gage network.

# Requirements
- Because the original data has a file for each 2-minute timestep, adding up to millions of files, an HPC was needed to download and process this data

# Inputs
- [The MRMS precipitation rate product ](https://vlab.noaa.gov/web/wdtd/-/surface-precipitation-rate-spr-?selectedFolder=9234881) was drawn from two sources:
  - 2001-2011:
    - Dr. Jian Zhang And Dr. Jonathan Gourley. (2018). Multi-Radar Multi-Sensor Precipitation Reanalysis (Version 1.0). Open Commons Consortium Environmental Data Commons. https://doi.org/10.25638/EDC.PRECIP.0001
  - 2015-present: [Iowa Environmental Mesonet](https://mesonet.agron.iastate.edu/)
    - Can be accessed through links like this: https://mtarchive.geol.iastate.edu/YYYY/MM/DD/mrms/ncep/PrecipRate/ (e.g., https://mtarchive.geol.iastate.edu/2020/09/09/mrms/ncep/PrecipRate/)

# Outputs
- Dataset
    - CONUS
        - `mrms_for_rainyday`
            - main output; contains a daily netcdf with the full resolution MRMS data
        - `mrms_for_rainyday_hourly`
            - daily netcdfs aggregated to an hourly timestep
        - `mrms_for_rainyday_daily`
            - daily snetcdfs aggregated to a daily timestep
        - `mrms_for_rainyday_daily_consolidated`
            - yearly netcdfs aggregated to daily timestep (i.e., same as `mrms_for_rainyday_daily` but combined in yearly netcdfs)
        - `mrms_zipped`
            - containis yearly zip folders containing the main output, `mrms_for_rainyday`
    - Spanning user defined gage network
        - `mrms_for_rainyday_subset_norfolk_csvs`
            - contains yearly csvs with mrms data at each gage location
        - `mrms_for_rainyday_subset_norfolk_netcdfs`
            - contains yearly netcdfs with the full resolution dataset clipped to a rectangle bounded by the maximum extend of the gage network plus a buffer of 5 gridcells in each direction
- Plots
    - Annual rainfall totals and annual rainfall anomolies (difference from mean)

# Plots
Note: for both plots, `robust=True`, meaning the 2nd and 98th percentiles were used for determining the color bar.

## Annual Rainfall Totals
![annual totals](https://github.com/lassiterdc/highres-radar-rainfall-processing/blob/main/plots/h_annual_statistics/all_years.png?raw=true)

## Annual Rainfall Anomolies (difference from the mean)
![annual anaomolies](https://github.com/lassiterdc/highres-radar-rainfall-processing/blob/main/plots/h_annual_statistics/all_years_anomolies_rel_to_avg_excluding_2012-2014.png?raw=true)


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
