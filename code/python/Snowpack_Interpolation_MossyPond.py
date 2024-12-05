import os
import pandas as pd
import arcpy

def process_station_data_by_year(station_file, station_name, apr1_avg, lat, lon, years, base_directory, date_col='DATE TIME', value_col='VALUE'):
    """
    Processes SWE data for multiple years and a given station.
    
    Args:
        station_file (str): File name for the station data.
        station_name (str): Three-letter station name.
        apr1_avg (float): April 1st average value for the station.
        lat (float): Latitude of the station.
        lon (float): Longitude of the station.
        years (list): List of years to process.
        base_directory (str): Base directory containing the station files.
        date_col (str): Column name for the date/time data.
        value_col (str): Column name for the SWE values.
    
    Returns:
        pd.DataFrame: Combined DataFrame for all years with additional information.
    """
    file_path = os.path.join(base_directory, station_file)
    data = pd.read_csv(file_path)
    data[date_col] = pd.to_datetime(data[date_col])
    
    # Generate April 1st dates for filtering
    april1_dates = [pd.Timestamp(f"{year}-04-01") for year in years]
    
    data_filtered = data[data[date_col].isin(april1_dates)].copy()
    data_filtered['Apr1Avg'] = apr1_avg
    data_filtered['pctAvg'] = round((data_filtered[value_col] / apr1_avg) * 100, 2)
    data_filtered['LAT'] = lat
    data_filtered['LON'] = lon
    data_filtered['Station'] = station_name
    
    return data_filtered

# List of years to process
years = list(range(2012, 2019))

# Base directory for SWE files
base_directory = r"C:\Mossy Pond\SWE_MossyPond\SWE_Data"

# Process each station
stations = [
    ("RDM.csv", "RDM", 48.7, 39.343, -120.508),
    ("DNS.csv", "DNS", 33.1, 39.31, -120.338),
    ("FOD.csv", "FOD", 36.5, 39.36, -120.5),
    ("CSL.csv", "CSL", 34.2, 39.325, -120.367),
    ("MWL.csv", "MWL", 52.8, 39.417, -120.508),
    ("ENM.csv", "ENM", 43.9, 39.436, -120.525),
    ("ONN.csv", "ONN", 20.1, 39.275, -120.358),
    ("BOM.csv", "BOM", 20.9, 39.458, -120.6),
    ("FNP.csv", "FNP", 29.8, 39.47, -120.572),
    ("WBB.csv", "WBB", 32.0, 39.485, -120.425)
]

# Combine data for all stations and years
all_data_frames = []
for station_file, station_name, apr1_avg, lat, lon in stations:
    station_data = process_station_data_by_year(
        station_file, station_name, apr1_avg, lat, lon, years, base_directory
    )
    all_data_frames.append(station_data)

# Combine all stations into one DataFrame
SWE_all = pd.concat(all_data_frames, ignore_index=True)
"C:\\Mossy Pond\\SWE_MossyPond\\Outputs\\SWE_MP_2014_2018.csv"
# Save combined data to CSV
SWE_all.to_csv(r"C:\Mossy Pond\SWE_MossyPond\Outputs\SWE_MP_2012_2018.csv", index=False)

# Process each year for IDW analysis
for year in years:
    # Subset data for the year
    SWE_year = SWE_all[SWE_all['DATE TIME'] == f"{year}-04-01"]
    year_csv_path = rf"C:\Mossy Pond\SWE_MossyPond\Outputs\SWE_MP_{year}.csv"
    SWE_year.to_csv(year_csv_path, index=False)

    # Convert CSV to standalone table
    arcpy.conversion.TableToTable(year_csv_path, r"C:\Mossy Pond\SWE_MossyPond\SWE_MossyPond.gdb", f"snowpack_table_{year}")
    
    # Create point feature class
    standalone_table = rf"C:\Mossy Pond\SWE_MossyPond\SWE_MossyPond.gdb\snowpack_table_{year}"
    output_fc = rf"C:\Mossy Pond\SWE_MossyPond\SWE_MossyPond.gdb\Snow_Stations_{year}"
    arcpy.management.XYTableToPoint(standalone_table, output_fc, "LON", "LAT", "pctAvg")

    # Run IDW
    IDW_snowpack = rf"C:\Mossy Pond\SWE_MossyPond\SWE_MossyPond.gdb\IDW_snowpack_{year}"
    IDW_result = arcpy.sa.Idw(in_point_features=output_fc, z_field="pctAvg")
    IDW_result.save(IDW_snowpack)

import os
import arcpy
import csv
from arcpy.sa import ZonalStatisticsAsTable

# Paths and Parameters
years = list(range(2012, 2019))  # Years to process
gdb_path = r"C:\Mossy Pond\SWE_MossyPond\SWE_MossyPond.gdb"  # Geodatabase path
study_area_fc = os.path.join(gdb_path, "MPSA_boundary")  # Updated MPSA Feature Class
zonal_stats_folder = r"C:\Mossy Pond\SWE_MossyPond\ZonalStats"  # Zonal stats outputs
output_csv = r"C:\Mossy Pond\SWE_MossyPond\Outputs\Average_Pixel_Values.csv"  # Output CSV

# Ensure the zonal stats folder exists
os.makedirs(zonal_stats_folder, exist_ok=True)

# Dictionary to store average pixel values
average_pixel_values = {}

# Process each year
for year in years:
    idw_raster = os.path.join(gdb_path, f"IDW_snowpack_{year}")
    zonal_stats_table = os.path.join(zonal_stats_folder, f"ZonalStats_{year}.dbf")
    
    # Perform Zonal Statistics
    ZonalStatisticsAsTable(
        in_zone_data=study_area_fc,
        zone_field="OBJECTID",
        in_value_raster=idw_raster,
        out_table=zonal_stats_table,
        statistics_type="MEAN"
    )
    
    # Extract the mean value
    with arcpy.da.SearchCursor(zonal_stats_table, ["MEAN"]) as cursor:
        for row in cursor:
            average_pixel_values[year] = row[0]
            break

# Write results to a CSV file
with open(output_csv, mode="w", newline="") as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(["Year", "Average Pixel Value"])
    for year, avg_value in average_pixel_values.items():
        writer.writerow([year, avg_value])



import os
import arcpy

# Paths and Parameters
gdb_path = r"C:\Mossy Pond\SWE_MossyPond\SWE_MossyPond.gdb"  # Geodatabase path
study_area_fc = os.path.join(gdb_path, "MPSA_boundary")  # Study Area Polygon Feature Class
stations_fc = os.path.join(gdb_path, "Snow_Stations_2018")  # Snow Stations Feature Class
centroid_fc = os.path.join(gdb_path, "MPSA_Centroid")  # Output Centroid Feature Class
out_table = os.path.join(gdb_path, "Centroid_To_Stations_Distances")  # Output Table
output_csv = r"C:\Mossy Pond\SWE_MossyPond\Outputs\Centroid_To_Snow_Stations.csv"  # Output CSV File

# Calculate the centroid of the MPSA boundary
arcpy.management.FeatureToPoint(study_area_fc, centroid_fc, point_location="CENTROID")

# Use PointDistance to calculate distances
arcpy.analysis.PointDistance(
    in_features=centroid_fc,
    near_features=stations_fc,
    out_table=out_table
)

# Export the output table to a CSV
arcpy.conversion.TableToTable(out_table, os.path.dirname(output_csv), os.path.basename(output_csv))




