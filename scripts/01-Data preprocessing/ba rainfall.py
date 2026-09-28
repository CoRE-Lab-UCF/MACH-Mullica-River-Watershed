#Loading packages
import xarray as xr
import numpy as np
import os
import glob
import pandas as pd

#Watershed name
name= 'Mullica'

# Read into a NumPy array
lat_lon_pairs = np.loadtxt(f'F:/OneDrive - University of Central Florida/Documents/CONUS/NJ/{name}/{name}_aorc_indicies.csv', 
                  delimiter=",")


#Number of AORC points
num_sites = len(lat_lon_pairs)

# AORC files
base_path = "F:/OneDrive - University of Central Florida/Documents/CONUS/Data/AORC/"
file_pattern = os.path.join(base_path, "AORC_1km_NE_*.nc")
aorc_files = sorted(glob.glob(file_pattern))

# Calculate total number of days across all years first
total_hours = 0
for file in aorc_files:
    year = os.path.basename(file).split('_')[-1].split('.')[0]
    # Open the file for the current year
    aorc_df = xr.open_dataset(file)
    
    # Calculate the number of days in this year
    year_start = pd.Timestamp(str(aorc_df['time'].min().values))
    year_end = pd.Timestamp(f"{year}-12-31 23:59:59")
    year_hours = int((year_end - year_start).total_seconds() // 3600) + 1
    total_hours += year_hours
    aorc_df.close()  # Close to free memory
    
print(f"Total hours across all years: {total_hours}")

# Pre-allocate the full array with the correct size
array = np.empty((num_sites, total_hours), dtype=np.float32)
array.fill(np.nan)  # Initialize with NaN

# Process the files
current_hour_idx = 0
for file in aorc_files:
    year = os.path.basename(file).split('_')[-1].split('.')[0]
    print(f"Processing {year}...")
    
    # Open the file for the current year
    aorc_df = xr.open_dataset(file)
    
    # Calculate the number of hours in this year
    year_start = pd.Timestamp(str(aorc_df['time'].min().values))
    year_end = pd.Timestamp(f"{year}-12-31 23:59:59") 
    year_hours = int((year_end - year_start).total_seconds() // 3600) + 1
    
    # Process each site
    for i, (lat_idx, lon_idx) in enumerate(lat_lon_pairs):
        # Extract hourly data
        series = aorc_df['APCP_surface'][:, int(lat_idx), int(lon_idx)]
        
        # Verify the length matches our expected hours in year
        if len(series) != year_hours:
            print(f"Warning: Expected {year_hours} hours for {year}, but got {len(series)}")
        
        # Store in the appropriate slice of the array
        hour_slice = slice(current_hour_idx, current_hour_idx + len(series))
        array[i, hour_slice] = series.values
    
    # Update the hours index for the next file
    current_hour_idx += year_hours
    
    # Close the dataset to free memory
    aorc_df.close()

print("Processing complete!")
print(f"Final array shape: {array.shape}")

np.savetxt(f'F:/OneDrive - University of Central Florida/Documents/CONUS/NJ/{name}/{name}_aorc_ba.csv', array, delimiter=",", fmt='%.6f')

##Calculate the average

import numpy as np

rain = np.loadtxt(f'C:/Users/ro327497/OneDrive - University of Central Florida/Documents/CONUS/NJ/{name}/{name}_aorc_ba.csv',  delimiter=',')

rain_avg = np.mean(rain, axis=0)  # axis=0 means average over rows (down columns)

# Result shape
print(rain_avg.shape)

np.savetxt(f'C:/Users/ro327497/OneDrive - University of Central Florida/Documents/CONUS/NJ/{name}/{name}_aorc_ba_ts.csv', rain_avg, delimiter=",", fmt='%.6f')





