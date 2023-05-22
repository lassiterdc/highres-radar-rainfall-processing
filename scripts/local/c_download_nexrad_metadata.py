#%% import libraries and load directories
from urllib.request import urlopen
import json
import pandas as pd
from datetime import datetime
import geopandas as gpd
from geodatasets import get_path
import matplotlib.pyplot as plt

from __filepaths import return_c_filepaths

url_stations, dir_data, f_out_station_data_with_dates, f_stations, f_nexrad_station_shp = return_c_filepaths()
#%% load station metadata so I can loop through all the stations to get the location and the start and end dates
df_stations = pd.read_csv(f_stations)

lst_sta_ids = []
lst_lats = []
lst_lons = []
lst_begin_dates = []
lst_end_dates = []
lst_still_collecting = []

for sta_id in df_stations["NCDCID"]:
    url = url_stations.format(sta_id)
    response = urlopen(url)
    dic_json = json.loads(response.read())

    lst_sta_ids.append(sta_id)
    lst_lats.append(dic_json["stationCollection"]['stations'][0]['header']['latitude_dec'])
    lst_lons.append(dic_json["stationCollection"]['stations'][0]['header']['longitude_dec'])
    lst_begin_dates.append(dic_json["stationCollection"]['stations'][0]['header']['por']['beginDate'])
    end_date = dic_json["stationCollection"]['stations'][0]['header']['por']['endDate']
    still_collecting = False
    if end_date.lower() == 'present':
        still_collecting = True
        end_date = str(datetime.now().date())
    lst_still_collecting.append(still_collecting)
    lst_end_dates.append(end_date)
#%% combine results into dataframe
df_station_data = pd.DataFrame(dict(station_id = lst_sta_ids,
                                    latitude = lst_lats,
                                    longitude = lst_lons,
                                    start_date = lst_begin_dates,
                                    end_date = lst_end_dates,
                                    still_collecting = lst_still_collecting))

df_station_data["start_date"] = pd.to_datetime(df_station_data.start_date)
df_station_data["end_date"] = pd.to_datetime(df_station_data.end_date)
df_station_data["record_lgnth_days"] = (df_station_data["end_date"] - df_station_data["start_date"]).dt.days
df_station_data["record_lgnth_years"] = df_station_data["record_lgnth_days"] / 365.25

df_station_data["start_date"] = (df_station_data.start_date).astype(str)
df_station_data["end_date"] = (df_station_data.end_date).astype(str)


df_station_data.to_csv(f_out_station_data_with_dates)

#%% create geopandas database
gdf = gpd.GeoDataFrame(
    df_station_data, geometry=gpd.points_from_xy(df_station_data.longitude, df_station_data.latitude), crs="EPSG:4326")

#%% test plotting
# world = gpd.read_file(get_path("naturalearth.land"))

# # We restrict to South America.
# ax = world.plot(color="white", edgecolor="black")

# # We can now plot our ``GeoDataFrame``.
# gdf.plot(ax=ax, color="red")

# plt.show()

#%% exporting as shapefile
gdf.to_file(f_nexrad_station_shp)