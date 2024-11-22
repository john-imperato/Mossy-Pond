import os
import pandas as pd

def process_station_data(station_file, station_name, apr1_avg, lat, lon, date_col='DATE TIME', value_col='VALUE'):
    """
    Processes snow water equivalent (SWE) data for CA snow survey stations.

    Args:
        station_file (str): File name for the station data.
        station_name (str): Three-letter station name.
        apr1_avg (float): April 1st average value for the station.
        lat (float): Latitude of the station.
        lon (float): Longitude of the station.
        date_col (str): Column name for the date/time data.
        value_col (str): Column name for the SWE values.

    Returns:
        pd.DataFrame: Processed DataFrame with added columns for April 1st average, pctAvg, and location.
    """
    # Define the base directory where your data is stored
    base_directory = r"C:\Modeling Geographical Objects\MossyPond\Data\SWE"
    
    # Join the base directory with the station file name to create the full file path
    file_path = os.path.join(base_directory, station_file)
    
    # Load the data and keep necessary columns
    data = pd.read_csv(file_path)
    
    # Filter only the relevant columns
    data = data[[date_col, value_col]]
    
    # Convert the date column to a datetime object
    data[date_col] = pd.to_datetime(data[date_col])
    
    # Define April 1st dates for filtering
    april1_dates = pd.to_datetime(["2012-04-01", "2013-04-01", "2014-04-01", "2015-04-01",
                                   "2016-04-01", "2017-04-01", "2018-04-01"])
    
    # Filter for April 1st observations
    data_filtered = data.loc[data[date_col].isin(april1_dates)].copy()
    
    # Add April 1st average column
    data_filtered.loc[:, 'Apr1Avg'] = apr1_avg
    
    # Calculate pctAvg (percentage of the April 1st average) and round to 2 decimals
    data_filtered.loc[:, 'pctAvg'] = round((data_filtered[value_col] / apr1_avg) * 100, 2)
    
    # Add location information
    data_filtered.loc[:, 'LAT'] = lat
    data_filtered.loc[:, 'LON'] = lon
    
    # Add station name column
    data_filtered.loc[:, 'Station'] = station_name
    
    return data_filtered

# RDM: Red Mountain
rdm_data = process_station_data("RDM.csv", "RDM", 48.7, 39.343, -120.508)

# DNS: Donner Summit
dns_data = process_station_data("DNS.csv", "DNS", 33.1, 39.31, -120.338)

# FOD: Fordyce Lake
fod_data = process_station_data("FOD.csv", "FOD", 36.5, 39.36, -120.5)

# CSL: Central Sierra Snow Lab
csl_data = process_station_data("CSL.csv", "CSL", 34.2, 39.325, -120.367)

# MWL: Meadow Lake
mwl_data = process_station_data("MWL.csv", "MWL", 52.8, 39.417, -120.508)

# ENM: English Mountain
enm_data = process_station_data("ENM.csv", "ENM", 43.9, 39.436, -120.525)

# ONN: Onion Creek
onn_data = process_station_data("ONN.csv", "ONN", 20.1, 39.275, -120.358)

# BOM: Bowman Lake
bom_data = process_station_data("BOM.csv", "BOM", 20.9, 39.458, -120.6)

# FNP: Findley Peak
fnp_data = process_station_data("FNP.csv", "FNP", 29.8, 39.47, -120.572)

# WBB: Webber Lake
wbb_data = process_station_data("WBB.csv", "WBB", 32.0, 39.485, -120.425)

# Combine data frames
data_frames = [rdm_data, dns_data, fod_data, csl_data, mwl_data, enm_data, onn_data, bom_data, fnp_data, wbb_data]
SWE_all = pd.concat(data_frames, ignore_index=True)

# save to .csv file
SWE_all.to_csv(r'C:\Modeling Geographical Objects\MossyPond\Data\SWE\SWE_MP.csv', index=False)

# subset data frame for 2017 winter
SWE_2017 = SWE_all[SWE_all['DATE TIME'] == "2017-04-01"]

# save subset to .csv file
SWE_2017.to_csv(r'C:\Modeling Geographical Objects\MossyPond\Data\SWE\SWE_MP_2017.csv', index=False)

SWE_2017

# merge all data frames
SWE_AllStations_ByYear = rdm_data.merge(dns_data, on='DATE TIME', how='outer', suffixes=('_RDM', '_DNS')) \
                      .merge(fod_data, on='DATE TIME', how='outer', suffixes=('', '_FOD')) \
                      .merge(csl_data, on='DATE TIME', how='outer', suffixes=('', '_CSL')) \
                      .merge(mwl_data, on='DATE TIME', how='outer', suffixes=('', '_MWL')) \
                      .merge(enm_data, on='DATE TIME', how='outer', suffixes=('', '_ENM')) \
                      .merge(onn_data, on='DATE TIME', how='outer', suffixes=('', '_ONN')) \
                      .merge(bom_data, on='DATE TIME', how='outer', suffixes=('', '_BOM')) \
                      .merge(fnp_data, on='DATE TIME', how='outer', suffixes=('', '_FNP')) \
                      .merge(wbb_data, on='DATE TIME', how='outer', suffixes=('', '_WBB'))




# define paths
csv_path = r"C:\Modeling Geographical Objects\MossyPond\Data\SWE\SWE_MP_2017.csv" # csv
gdb_path = r"C:\Modeling Geographical Objects\MossyPond\GDB_MP.gdb" # geodatabase
standalone_table = r"C:\Modeling Geographical Objects\MossyPond\GDB_MP.gdb\snowpack_table"  # standalone table
output_fc = r"C:\Modeling Geographical Objects\MossyPond\GDB_MP.gdb\Snow_Stations" # feature class

# Create standalone table
arcpy.conversion.TableToTable(csv_path, gdb_path, "snowpack_table")

# create point feature class
arcpy.management.XYTableToPoint(standalone_table, output_fc, "LON", "LAT", "pctAvg")

# IDW

# set environment extent
left = -120.687559785412
right = -120.280688043732
bottom = 39.2392347691052
top = 39.5243532244492
arcpy.env.extent = arcpy.Extent(left, bottom, right, top)

# input parameters
z_field = "pctAvg"
IDW_snowpack = r"C:\Modeling Geographical Objects\MossyPond\GDB_MP.gdb\IDW_snowpack"

# run IDW
IDW_sno = arcpy.sa.Idw(
    in_point_features=output_fc,
    z_field=z_field,
)

# save output raster
IDW_sno.save(r"C:\Modeling Geographical Objects\MossyPond\GDB_MP.gdb\IDW_snowpack_2017")
