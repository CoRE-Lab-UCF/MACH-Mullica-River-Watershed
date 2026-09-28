import xarray as xr
import numpy as np
import os
import glob
import pandas as pd

# Precompute lat-lon pairs as tuples
lat_lon_pairs = list(zip(lat_indices, lon_indices))
num_sites = len(lat_lon_pairs)

# AORC files
base_path = "F:/OneDrive - University of Central Florida/Documents/CONUS/Data/AORC/"
file_pattern = os.path.join(base_path, "AORC_1km_NE_*.nc")
aorc_files = sorted(glob.glob(file_pattern))

# Calculate total number of days across all years first
total_days = 0
for file in aorc_files:
    year = os.path.basename(file).split('_')[-1].split('.')[0]
    # Open the file for the current year
    aorc_df = xr.open_dataset(file)
    
    # Calculate the number of days in this year
    year_start = pd.Timestamp(str(aorc_df['time'].min().values))
    year_end = pd.Timestamp(f"{year}-12-31")
    year_days = (year_end - year_start).days + 1
    total_days += year_days

print(f"Total days across all years: {total_days}")

# Pre-allocate the full array with the correct size
daily_array = np.empty((num_sites, total_days), dtype=np.float32)
daily_array.fill(np.nan)  # Initialize with NaN

# Process the files
current_day_idx = 0
for file in aorc_files:
    year = os.path.basename(file).split('_')[-1].split('.')[0]
    print(f"Processing {year}...")
    
    # Open the file for the current year
    aorc_df = xr.open_dataset(file)
    
    # Calculate the number of days in this year
    year_start = pd.Timestamp(str(aorc_df['time'].min().values))
    year_end = pd.Timestamp(f"{year}-12-31")
    year_days = (year_end - year_start).days + 1
    
    # Process each site
    for i, (lat_idx, lon_idx) in enumerate(lat_lon_pairs):
        # Extract and resample hourly data to daily
        daily_series = aorc_df['APCP_surface'][:, lat_idx, lon_idx].resample(time='1D').sum()
        
        # Verify the length matches our expected days in year
        if len(daily_series) != year_days:
            print(f"Warning: Expected {year_days} days for {year}, but got {len(daily_series)}")
        
        # Store in the appropriate slice of the array
        day_slice = slice(current_day_idx, current_day_idx + len(daily_series))
        daily_array[i, day_slice] = daily_series.values
    
    # Update the day index for the next file
    current_day_idx += year_days
    
    # Close the dataset to free memory
    aorc_df.close()

print("Processing complete!")
print(f"Final array shape: {daily_array.shape}")

np.savetxt(f'F:/OneDrive - University of Central Florida/Documents/CONUS/NJ/{name}/{name}_aorc_.csv', daily_array, delimiter=",", fmt='%.6f')
df = pd.DataFrame(lat_lon_pairs)
df.to_csv(f'F:/OneDrive - University of Central Florida/Documents/CONUS/NJ/{name}/{name}_aorc_indicies.csv', index=False, header=False)
