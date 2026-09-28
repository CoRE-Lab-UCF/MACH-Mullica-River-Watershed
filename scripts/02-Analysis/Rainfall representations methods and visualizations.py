# Import packages
import pandas as pd
import matplotlib.pyplot as plt
import geopandas as gpd
from shapely.geometry import Point
import numpy as np
import cartopy.io.img_tiles as cimgt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import xarray as xr
import os
import glob
from pathlib import Path

# Basin name
name = 'Mullica'

# Read in cora centroids
cora = pd.read_csv(f'C:/Users/ro327497/OneDrive - University of Central Florida/Documents/CONUS/NJ/{name}/{name}_cora_centroid_lat_lon.csv')

# AORC points in the watershed
aorc = pd.read_csv(f'C:/Users/ro327497/OneDrive - University of Central Florida/Documents/CONUS/NJ/{name}/{name}_aorc_subset.csv')

# Indicies of AORC points in the watershed
aorc_indicies = pd.read_csv(f'C:/Users/ro327497/OneDrive - University of Central Florida/Documents/CONUS/NJ/{name}/{name}_aorc_indicies.csv', header = None, names=['lat_idx', 'lon_idx'])

#Distances to test in km
distances = [5,10,15,20,25,30]

##Identify AORC points within these disances of 

# Initialize empty dataframe to store all results
all_results = pd.DataFrame()

for dist in distances:
 
 # Convert CORA to GeoDataFrame
 cora_gdf = gpd.GeoDataFrame(cora, geometry=gpd.points_from_xy(cora.lon, cora.lat),crs="EPSG:4326")

 # Convert AORC to GeoDataFrame
 aorc_gdf = gpd.GeoDataFrame(aorc, geometry=gpd.points_from_xy(aorc.lon, aorc.lat),crs="EPSG:4326")

 # Convert to projected CRS (UTM Zone 18N for New Jersey)
 cora_gdf_proj = cora_gdf.to_crs(epsg=32618)  # UTM 18N in meters
 aorc_gdf_proj = aorc_gdf.to_crs(epsg=32618)

 # Create buffers of 5 km (5000 meters)
 cora_gdf_proj['buffer'] = cora_gdf_proj.geometry.buffer(dist * 1000)

 # Clean up any existing index columns
 if 'index_right' in aorc_gdf_proj.columns:
     aorc_gdf_proj = aorc_gdf_proj.drop(columns=['index_right'])
 if 'index_right' in cora_gdf_proj.columns:
     cora_gdf_proj = cora_gdf_proj.drop(columns=['index_right'])

 # Spatial join AORC points within buffer polygons
 joined = gpd.sjoin(
    aorc_gdf_proj, 
    cora_gdf_proj[['buffer', 'cluster_id']].set_geometry('buffer'), 
    predicate='within',
    how='inner'
 )

 print(f"Number of AORC points matched: {len(joined)}")
 print(joined.head())


 # Convert back to WGS84 for plotting
 joined = joined.to_crs(epsg=4326)
 print(joined)
 # Get aorc_df indicies of lat and lon points   
 joined['lat_idx'] = aorc_indicies.loc[joined.index, 'lat_idx'].values
 joined['lon_idx'] = aorc_indicies.loc[joined.index, 'lon_idx'].values
 
 joined = joined[['lon','lat','lon_idx','lat_idx','geometry','DEPWMAS_ID','WMA','WREGION','cluster_id']]
 joined = joined.sort_values('cluster_id').reset_index(drop=True)
 joined.to_csv(f'C:/Users/ro327497/OneDrive - University of Central Florida/Documents/CONUS/NJ/{name}/Wahl approach/{name}_aorc_centroids_within_{dist}_km.csv', index=False, header=True)

 # Define labels
 labels = ['A', 'B', 'C', 'D', 'E']

 # Define terrain
 stamen_terrain = cimgt.StadiaMapsTiles(
     apikey='1dabf9e0-8140-4bfc-9327-2c1ef91e4557', 
     style='stamen_terrain_background', 
     resolution='@2x'
 ) 

 # Create the figure
 fig = plt.figure(figsize=(15, 3.5))

 # Create subplot grid
 gs = fig.add_gridspec(1, 5, wspace=0.01)

 # Loop through 5 CORA centroids
 for i in range(5):
    ax = fig.add_subplot(gs[0, i], projection=stamen_terrain.crs)
    
    # Set map extent (adjust as needed)
    ax.set_extent([-75, -74, 39, 40], crs=ccrs.PlateCarree())
    
    # Add terrain and features
    ax.add_image(stamen_terrain, 8)
    ax.add_feature(cfeature.COASTLINE)
    ax.add_feature(cfeature.STATES)
    ax.add_feature(cfeature.RIVERS)
    ax.add_feature(cfeature.LAKES)
    ax.add_feature(cfeature.OCEAN)
    ax.add_feature(cfeature.LAND)
    
    # Get subset for this cluster
    subset = joined[joined["cluster_id"] == i]
    
    print(f"Cluster {i+1}: {len(subset)} points")
    
    # Plot AORC points in red (small dots)
    if len(subset) > 0:
        ax.scatter(subset.lon, subset.lat, s=0.1, c="Red", transform=ccrs.PlateCarree())
    
    # Plot the CORA centroid point in black
    point = cora.iloc[i]
    ax.scatter(point.lon, point.lat, color='black', s=20, 
               edgecolors='black', linewidths=1, transform=ccrs.PlateCarree())
    
    # Add title with label
    ax.set_title(f'CORA centroid {labels[i]}', fontsize=14)

    if i == 0:
        ax.text(-0.175, 0.5, f'{dist} km', transform=ax.transAxes,
        va='center', ha='center', fontsize=14)
        
 # Adjust layout and display
 plt.tight_layout(rect=[0, 0, 0.9, 1])

 # Save figure
 plt.savefig(f'C:/Users/ro327497/OneDrive - University of Central Florida/Documents/CONUS/NJ/{name}/Wahl approach/anv_{name}_within_{dist}_km.png', dpi=300, bbox_inches='tight')
 plt.show()


joined.index


# Get average hourly series for each cora point and distance

# Define labels
labels = ['A', 'B', 'C', 'D', 'E']

# distances
distances = [10, 15, 20, 25, 30]

for dist in distances:
 #Loop through distances and cora points
 for cluster_id  in range(5):
    
  # Read in lat/lon of aorc points within dist
  lat_lon_pairs = pd.read_csv(f'C:/Users/ro327497/OneDrive - University of Central Florida/Documents/CONUS/NJ/{name}/Wahl approach/{name}_aorc_centroids_within_{dist}_km.csv')
  lat_lon_pairs = lat_lon_pairs[lat_lon_pairs['cluster_id'] == cluster_id ]
  lat_lon_pairs = lat_lon_pairs[['lat_idx','lon_idx']]
  print(lat_lon_pairs.head())
  #Number of AORC points
  num_sites = len(lat_lon_pairs)

  # AORC files
  base_path = "C:/Users/ro327497/OneDrive - University of Central Florida/Documents/CONUS/Data/AORC/"
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
    for i, (lat_idx, lon_idx) in enumerate(lat_lon_pairs.itertuples(index=False)):
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

  # Average rainfall by hour
  rain  = array.mean(axis=0)  # axis=0 means average over rows (down columns)
 
  # Write to CSV file
  pd.Series(rain).to_csv(f'C:/Users/ro327497/OneDrive - University of Central Florida/Documents/CONUS/NJ/{name}/Wahl approach/{name}_cora_centroid_{labels[cluster_id]}_{dist}_km_aorc_ba.csv', index=False,  header = False)

# Check lat/lons work
aorc_df = xr.open_dataset(file)
aorc_df['longitude'][1880]
aorc_df['latitude'][775]


# Find the controid of each set of aorc points
# Read in lat/lon of aorc points within dist

for dist in distances:

    aorc_points = pd.read_csv(f'C:/Users/ro327497/OneDrive - University of Central Florida/Documents/CONUS/NJ/{name}/Wahl approach/{name}_aorc_centroids_within_{dist}_km.csv')

    aorc_centroids = []

    for cluster_id in range(5):

        # Select one cluster
        cluster = aorc_points[aorc_points["cluster_id"] == cluster_id]

        # Create GeoDataFrame
        gdf = gpd.GeoDataFrame(cluster, geometry=gpd.points_from_xy(cluster["lon"], cluster["lat"]), crs="EPSG:4326")

        # Project to UTM
        gdf = gdf.to_crs(32618)

        # Compute centroid
        centroid = gdf.unary_union.centroid

        # Convert back to WGS84
        centroid = gpd.GeoSeries([centroid], crs=32618).to_crs(4326)

        # Save coordinates
        aorc_centroids.append({"cluster_id": cluster_id, "lon": centroid.x.iloc[0], "lat": centroid.y.iloc[0]})

    # Convert to DataFrame and save
    aorc_centroids = pd.DataFrame(aorc_centroids)

    aorc_centroids.to_csv(f'C:/Users/ro327497/OneDrive - University of Central Florida/Documents/CONUS/NJ/{name}/Wahl approach/{name}_centroid_aorc_within_{dist}_km.csv', index=False)



# Centroid for overall basin average series
aorc_points = pd.read_csv(f'C:/Users/ro327497/OneDrive - University of Central Florida/Documents/CONUS/NJ/{name}/{name}_aorc_subset.csv')

# Create GeoDataFrame
gdf = gpd.GeoDataFrame(aorc_points.index_right, geometry=gpd.points_from_xy(aorc_points["lon"], aorc_points["lat"]), crs="EPSG:4326")

# Project to UTM
gdf = gdf.to_crs(32618)

# Compute centroid
centroid = gdf.unary_union.centroid

# Convert back to WGS84
centroid = gpd.GeoSeries([centroid], crs=32618).to_crs(4326)

# Convert coords to DataFrame and save
aorc_centroids = pd.DataFrame({"lon": [centroid.x.iloc[0]], "lat": [centroid.y.iloc[0]]})

aorc_centroids.to_csv(f'C:/Users/ro327497/OneDrive - University of Central Florida/Documents/CONUS/NJ/{name}/Wahl approach/{name}_subset_centroid_aorc.csv', index=False)




# Showing points within given distance

# Read in an AORC file
aorc_df = xr.open_dataset('C:/Users/ro327497/OneDrive - University of Central Florida/Documents/CONUS/Data/AORC\\AORC_1km_NE_1993.nc')

# Read in aorc points located within distances of cora points
dist = 25
joined = pd.read_csv(f'C:/Users/ro327497/OneDrive - University of Central Florida/Documents/CONUS/NJ/{name}/Wahl approach/{name}_aorc_centroids_within_{dist}_km.csv')

# Define labels
labels = ['A', 'B', 'C', 'D', 'E']

# Define terrain
stamen_terrain = cimgt.StadiaMapsTiles(
    apikey='1dabf9e0-8140-4bfc-9327-2c1ef91e4557', 
    style='stamen_terrain_background', 
    resolution='@2x'
) 

# Create the figure
fig = plt.figure(figsize=(15, 3.5))

# Create subplot grid
gs = fig.add_gridspec(1, 5, wspace=0.01)

# Loop through 5 CORA centroids
for i in range(5):
    ax = fig.add_subplot(gs[0, i], projection=stamen_terrain.crs)
    
    # Set map extent (adjust as needed)
    ax.set_extent([-75, -74, 39, 40], crs=ccrs.PlateCarree())
    
    # Add terrain and features
    ax.add_image(stamen_terrain, 8)
    ax.add_feature(cfeature.COASTLINE)
    ax.add_feature(cfeature.STATES)
    ax.add_feature(cfeature.RIVERS)
    ax.add_feature(cfeature.LAKES)
    ax.add_feature(cfeature.OCEAN)
    ax.add_feature(cfeature.LAND)
    
    # Get subset for this cluster
    subset = joined[joined["cluster_id"] == i]
    subset = subset[['lon_idx', 'lat_idx']]
    
    print(f"Cluster {i+1}: {len(subset)} points")
    
    # Plot AORC points in red (small dots)
    if len(subset) > 0:
       lon_vals = aorc_df['longitude'].values  # shape (ny, nx)
       lat_vals = aorc_df['latitude'].values
       ax.scatter(lon_vals[subset.lon_idx.values.astype(int)],
                  lat_vals[subset.lat_idx.values.astype(int)],
                  s=0.1, c="Red", transform=ccrs.PlateCarree())
    
    # Plot the CORA centroid point in black
    point = cora.iloc[i]
    ax.scatter(point.lon, point.lat, color='black', s=20, 
               edgecolors='black', linewidths=1, transform=ccrs.PlateCarree())
    
    # Add title with label
    ax.set_title(f'CORA centroid {labels[i]}', fontsize=12)

# Adjust layout and display
plt.tight_layout(rect=[0, 0, 0.9, 1])

# Close the dataset to free memory
aorc_df.close()

# plot figure
plt.show()


# Read in results
dfs = []
rain_dur = 24
distances = [5, 10, 15, 20, 25, 30]

for dist in distances:
 for cora in labels:
  df =  pd.read_csv(f'C:/Users/ro327497/OneDrive - University of Central Florida/Documents/CONUS/NJ/{name}/Wahl Approach/{name}_anv_spatial_{rain_dur}_hr_acc_time_10_yr_r_ntr_cora_{cora}_dist_{dist}.csv',
                    index_col=False)
  dfs.append(df)

# Concatenate by columns (axis=1)
res_df = pd.concat(dfs, axis=0, ignore_index=True)

##Plot to represent method 2
columns_to_analyze = ['cor_r_tc', 'cor_r_non_tc',
                       'cor_ntr_tc', 'cor_ntr_non_tc', 
                       'cor_pval_r_tc', 'cor_pval_r_non_tc', 'cor_pval_ntr_tc',
                       'cor_pval_ntr_non_tc']

results = {}

for col in columns_to_analyze:
 cora_mean = []
 cora_sd = []
 cora_p5 = []   # 5th percentile
 cora_p50 = []  # Median
 cora_p95 = []  # 95th percentile

 for i in range(5):
  # Get data for this CORA cluster
  cora_data = res_df.loc[res_df['cora'] == i+1, col]

  # Calculate statistics
  ave = cora_data.mean()
  std_dev = np.std(cora_data, ddof=1)
  p5 = np.percentile(cora_data, 5)
  p50 = np.percentile(cora_data, 50)
  p95 = np.percentile(cora_data, 95)

  # Append to lists
  cora_mean.append(ave)
  cora_sd.append(std_dev)
  cora_p5.append(p5)
  cora_p50.append(p50)
  cora_p95.append(p95)

 results[col] = {'mean': cora_mean, 'std': cora_sd, 'p5': cora_p5, 'p50': cora_p50, 'p95': cora_p95}
    
    
# Read in results
dfs = []
rain_dur = 24
distances = [5, 10, 15, 20, 25, 30]

for dist in distances:
 for cora in labels:
  df =  pd.read_csv(f'C:/Users/ro327497/OneDrive - University of Central Florida/Documents/CONUS/NJ/{name}/Wahl Approach/{name}_anv_spatial_{rain_dur}_hr_acc_time_10_yr_r_ntr_cora_{cora}_dist_{dist}.csv',
                    index_col=False)
  dfs.append(df)

# Concatenate by columns (axis=1)
res_df = pd.concat(dfs, axis=0, ignore_index=True)

columns_to_plot = ['cor_r_tc', 'cor_r_non_tc','cor_ntr_tc', 'cor_ntr_non_tc']
p_val_cols = ['cor_pval_r_tc', 'cor_pval_r_non_tc', 'cor_pval_ntr_tc','cor_pval_ntr_non_tc']
col_to_pval = {
    'cor_r_tc': 'cor_pval_r_tc',
    'cor_r_non_tc': 'cor_pval_r_non_tc',
    'cor_ntr_tc': 'cor_pval_ntr_tc',
    'cor_ntr_non_tc': 'cor_pval_ntr_non_tc',
}

for col in columns_to_plot:
 #Create the figure
 fig = plt.figure(figsize=(15, 3.5))
 
 # Create subplot grid for manual assignment
 gs = fig.add_gridspec(1, 5, wspace=0.01)

 # Create each subplot individually with the correct projection
 axes = []
 for i in range(5):  # Create 5 subplots (3×2 grid minus 1)
  row = 0 #i // 2
  grid_col  = i #% 2
 
  #Create the plot
  ax = fig.add_subplot(gs[row, grid_col ])
 
  subset =  res_df[res_df['cora']==(i+1)]
  significant = subset[col_to_pval[col]] < 0.05
  
  # Median line
  #ax.plot(distances, subset[col], marker='o',  label=f'cora {i}')
  # Line itself (no markers, or markers styled separately below) 
  ax.plot(distances, subset[col], color='C0', zorder=1)

  # Significant points: filled
  ax.scatter(np.array(distances)[significant], subset[col][significant],
           marker='o', facecolors='C0', edgecolors='C0', zorder=2)

  # Non-significant points: hollow/open
  ax.scatter(np.array(distances)[~significant], subset[col][~significant],
           marker='o', facecolors='none', edgecolors='C0', zorder=2)

  # X-axis settings: fixed limits and ticks
  ax.set_xlim(0, 35)
  ax.set_xticks([1, 5, 10, 15, 20, 25, 30])  # ticks every 10 units
  ax.set_xlabel('Buffer Distance (km)')
   
  ax.set_ylim(-0.5, 0.5)
 
  # Y-axis only for first subplot
  if i == 0:
     ax.set_ylabel('Correlation')
  else:
     ax.set_yticks([])
        
  # Horizontal lines for mean ± 1 std
  #mean_val = cora_ave[i] mean_val + std_dev
  #std_dev = cora_sd[i]
  ax.axhline(results[col]['p50'][i], color='black', linestyle='--', linewidth=1, label='Mean')
  ax.axhline(results[col]['p5'][i], color='black', linestyle=':', linewidth=1, label='+1 SD')
  ax.axhline(results[col]['p95'][i], color='black', linestyle=':', linewidth=1, label='-1 SD')
 
  # Title
  ax.set_title(f'CORA centroid {labels[i]}', fontsize=12)

 # Display the plot
 plt.tight_layout()

 #save plot
 plt.savefig(f'C:/Users/ro327497/OneDrive - University of Central Florida/Documents/CONUS/NJ/{name}/Wahl approach/{name}_distance_sensitivity_{col}.png',bbox_inches='tight',pad_inches=0.1)

 plt.show()






#Plotting data by distance (Wahl method)

#Finding the confidence intervals based on all the data
# Read csvs and concatenate by columns
dfs = []
files = [500, 1000, 1500, 2000, 2300]
for i in files:
   Name =  f'{name}_anv_spatial_24_hr_acc_time_10_yr_r_ntr_{i}'
   df = pd.read_csv(
        f'C:/Users/ro327497/OneDrive - University of Central Florida/Documents/Paper/prav/results/con_ntr/data/{Name}.csv',
        index_col=False#, 
        #usecols=[0, 1, 2, 3, 4, 5, 18, 19, 20, 21, 31]
    )
   dfs.append(df)

# Concatenate by columns (axis=1)
data = pd.concat(dfs, axis=0, ignore_index=True)


##Plot to represent method 2
cora_ave = []
cora_sd = []
cora_p5 = []   # 5th percentile
cora_p50 = []   # Median
cora_p95 = []  # 95th percentile

for i in range(5):
    # Get data for this CORA cluster
    cora_data =np.minimum(data.loc[data['cora'] == i+1, 'rp_joint'], 100)
    
    # Calculate statistics
    ave = cora_data.mean()
    std_dev = np.std(cora_data, ddof=1)
    p5 = np.percentile(cora_data, 5)
    p50 = np.percentile(cora_data, 50)
    p95 = np.percentile(cora_data, 95)
    
    # Append to lists
    cora_ave.append(ave)
    cora_sd.append(std_dev)
    cora_p5.append(p5)
    cora_p50.append(p50)
    cora_p95.append(p95)
    
#Range of 95% confidence intervals
# Range of 5–95 percentile interval
cora_range = np.array(cora_p95) - np.array(cora_p5)

print(cora_range)

# Read in results
dfs = []
rain_dur = 24
distances = [5, 10, 15, 20, 25, 30]

for dist in distances:
 for cora in labels:
  df =  pd.read_csv(f'C:/Users/ro327497/OneDrive - University of Central Florida/Documents/CONUS/NJ/{name}/Wahl Approach/{name}_anv_spatial_{rain_dur}_hr_acc_time_10_yr_r_ntr_cora_{cora}_dist_{dist}.csv',
                    index_col=False)
  dfs.append(df)

# Concatenate by columns (axis=1)
res_df = pd.concat(dfs, axis=0, ignore_index=True)

#Create the figure
fig = plt.figure(figsize=(15, 3.5))
 
# Create subplot grid for manual assignment
gs = fig.add_gridspec(1, 5, wspace=0.01)

# Create each subplot individually with the correct projection
axes = []
for i in range(5):  # Create 5 subplots (3×2 grid minus 1)
 row = 0 #i // 2
 col = i #% 2
 
 #Create the plot
 ax = fig.add_subplot(gs[row, col])
 
 subset =  res_df[res_df['cora']==(i+1)]
 
 # Median line
 ax.plot(distances, np.minimum(subset['rp_joint'],100), marker='o',  label=f'cora {i}')
 
 # X-axis settings: fixed limits and ticks
 ax.set_xlim(0, 35)
 ax.set_xticks([1, 5, 10, 15, 20, 25, 30])  # ticks every 10 units
 ax.set_xlabel('Buffer Distance (km)')
   
 ax.set_ylim(20, 105)
 
 # Y-axis only for first subplot
 if i == 0:
     ax.set_ylabel('Joint return period (years)')
 else:
     ax.set_yticks([])
        
 # Horizontal lines for mean ± 1 std
 #mean_val = cora_ave[i] mean_val + std_dev
 #std_dev = cora_sd[i]
 ax.axhline(cora_p50[i], color='black', linestyle='--', linewidth=1, label='Mean')
 ax.axhline(cora_p95[i], color='black', linestyle=':', linewidth=1, label='+1 SD')
 ax.axhline(cora_p5[i], color='black', linestyle=':', linewidth=1, label='-1 SD')
 
 # Title
 ax.set_title(f'CORA centroid {labels[i]}', fontsize=12)

# Display the plot
plt.tight_layout()

#save plot
plt.savefig(f'C:/Users/ro327497/OneDrive - University of Central Florida/Documents/CONUS/NJ/{name}/Wahl approach/{name}_distance_sensitivity_10_yr.png',bbox_inches='tight',pad_inches=0.1)

plt.show()











#AORC points within 25km of CORA points

#Cora points
cora_clean = pd.read_csv(f'C:/Users/ro327497/OneDrive - University of Central Florida/Documents/CONUS/NJ/{name}/{name}_cora_clean.csv')
cora_centroids_gdf = pd.read_csv(f'C:/Users/ro327497/OneDrive - University of Central Florida/Documents/CONUS/NJ/{name}/{name}_cora_centroid_lat_lon.csv')

#AORC points in the watershed
aorc = pd.read_csv(f'C:/Users/ro327497/OneDrive - University of Central Florida/Documents/CONUS/NJ/{name}/{name}_aorc_subset.csv')

#Read in dataframe containing closest AORC points to CORA centroids
closest_points = pd.read_csv( f'C:/Users/ro327497/OneDrive - University of Central Florida/Documents/CONUS/NJ/{name}/{name}_closest_aorc_centroid_to_each_cora.csv', header=0)

#Define labels
labels = ['A','B','C','D','E']

# Read in results  of method 2 using 30km radius
dfs = []
rain_dur = 24
dist = 30

for cora in labels:
 df =  pd.read_csv(f'C:/Users/ro327497/OneDrive - University of Central Florida/Documents/CONUS/NJ/{name}/Wahl Approach/{name}_anv_spatial_{rain_dur}_hr_acc_time_10_yr_r_ntr_cora_{cora}_dist_{dist}.csv',
                   index_col=False)
 dfs.append(df)

# Concatenate by columns (axis=1)
res_df = pd.concat(dfs, axis=0, ignore_index=True)

#Read in basin average for each cora point
dfs = []
for cora in labels:
 df =  pd.read_csv(f'C:/Users/ro327497/OneDrive - University of Central Florida/Documents/CONUS/NJ/{name}/Wahl Approach/{name}_anv_spatial_24_hr_acc_time_10_yr_r_ntr_cora_{cora}_ba.csv',
                   index_col=False)
 dfs.append(df)

# Concatenate by columns (axis=1)
ba_ave = pd.concat(dfs, axis=0, ignore_index=True)


#Define terrain
stamen_terrain = cimgt.StadiaMapsTiles(apikey='1dabf9e0-8140-4bfc-9327-2c1ef91e4557', style='stamen_terrain_background', resolution='@2x')

#Create the figure
fig = plt.figure(figsize=(15, 3.5))
 
# Create subplot grid for manual assignment
gs = fig.add_gridspec(1, 5, wspace=0.01)

# Create each subplot individually with the correct projection
axes = []

#Loopthrouhg different number of events per year
for i in range(5):  # Create 5 subplots (3×2 grid minus 1)
 row = 0 #i // 2
 col = i #% 2
 ax = fig.add_subplot(gs[row, col], projection=stamen_terrain.crs)
 #ax = fig.add_axes(positions[i], projection=stamen_terrain.crs)
 axes.append(ax)

 cmap = plt.get_cmap('plasma_r')
 # set miami-dade county as map extent
 ax.set_extent([-75, -74, 39, 40], crs=ccrs.PlateCarree())
 # Add the Stamen data at zoom level 8.

 ax.add_image(stamen_terrain, 8)
 ax.add_feature(cfeature.COASTLINE)
 ax.add_feature(cfeature.STATES)
 ax.add_feature(cfeature.RIVERS)
 ax.add_feature(cfeature.LAKES)
 ax.add_feature(cfeature.OCEAN)
 ax.add_feature(cfeature.LAND)
 
 #Pull out lat lon of AORC index   
 idx = int(closest_points.aorc_index[i])
 lon = aorc.iloc[idx].lon
 lat = aorc.iloc[idx].lat
 
 # Use scatter to control marker size and transformation
 ax.scatter(lon, lat, s=10, c = "Red",  transform=ccrs.PlateCarree())

 #Add CORA point
 point = cora_centroids_gdf.iloc[i]
 ax.scatter(point.lon, point.lat, color='black', s=10, edgecolors='black', linewidths=1, transform=ccrs.PlateCarree())

 # Add title with label
 ax.set_title(f'CORA centroid {labels[i]}', fontsize=12)
 
# Display the plot
plt.tight_layout(rect=[0, 0, 0.9, 1])  # Adjust layout but leave space for colorbar

#save plot
plt.savefig(f'C:/Users/ro327497/OneDrive - University of Central Florida/Documents/CONUS/NJ/{name}/anv_{name}_cloest_aorc_points.png')

plt.show()

#Read in dataframe containing closest AORC points to CORA centroids
closest_points = pd.read_csv( f'C:/Users/ro327497/OneDrive - University of Central Florida/Documents/CONUS/NJ/{name}/{name}_closest_aorc_centroid_to_each_cora.csv', header=0)

data[(data['cora']==1) & (data['aorc']==int(closest_points.aorc_index[0]))] #45.241925
data[(data['cora']==2) & (data['aorc']==int(closest_points.aorc_index[1]))] #29.283941
data[(data['cora']==3) & (data['aorc']==int(closest_points.aorc_index[2]))] #37.313358
data[(data['cora']==4) & (data['aorc']==int(closest_points.aorc_index[3]))] #35.075455
data[(data['cora']==5) & (data['aorc']==int(closest_points.aorc_index[4]))] #29.419686

point_est = [45.241925, 29.283941, 37.313358, 35.075455, 29.419686]

# Assuming you have these lists already:
# cora_ave, cora_p5, cora_p95, and another_estimate (your cross markers)

# Create figure
fig, ax = plt.subplots(figsize=(10, 6))

# X-axis: CORA centroids
cora_ids = np.arange(1, 6)
labels = ['A', 'B', 'C', 'D', 'E']

# Plot shaded area between 5th and 95th percentiles
for i, cora_id in enumerate(cora_ids):
    ax.plot([cora_id, cora_id], [cora_p5[i], cora_p95[i]], 
            color='steelblue', linewidth=2, alpha=0.5, zorder=1)
    # Optional: Add caps at the ends
    ax.plot([cora_id-0.075, cora_id+0.075], [cora_p5[i], cora_p5[i]], 
            color='steelblue', linewidth=2, alpha=0.5)
    ax.plot([cora_id-0.075, cora_id+0.075], [cora_p95[i], cora_p95[i]], 
            color='steelblue', linewidth=2, alpha=0.5)
    ax.plot([cora_id-0.07, cora_id+0.07], [cora_p50[i], cora_p50[i]], 
            color='darkblue', linewidth=2, alpha=0.6)


# Add a single label for the legend — observed data (all CORA grid points)
ax.plot([], [], color='steelblue', linewidth=2, alpha=0.5, 
        label='Catchment estimates: 5th–95th percentile')
ax.plot([], [], color='darkblue', linewidth=2, alpha=0.6, 
        label='Catchment estimates: median')

# Plot mean as dots
ax.plot(cora_ids, cora_ave, 'o', markersize=12, 
      color='darkblue', label='Catchment estimates: mean', zorder=3)

# Composite record estimates with radius 30km (Method 2)
ax.plot(cora_ids, np.minimum(res_df['rp_joint'], 100), 'D', markersize=10,  
        color= 'Green', label='Method 2 - Local composite (30 km radius)', zorder=3)

# Plot ba mean as dots — Method 3 (basin average)
ax.plot(cora_ids, ba_ave['rp_joint'], '^', markersize=10, 
        color='#E69F00', label='Method 3 - Basin average', zorder=3)

# Plot point-to-point estimates as crosses — Method 1
ax.plot(cora_ids, point_est, 'x', markersize=12, 
        markeredgewidth=3, color='red', 
        label='Method 1 - Point-to-point', zorder=3)

handles, label = ax.get_legend_handles_labels()

desired_order = [
    'Catchment estimates: mean',
    'Catchment estimates: median',
    'Catchment estimates: 5th–95th percentile',
    'Method 1 - Point-to-point',
    'Method 2 - Local composite (30 km radius)',
    'Method 3 - Basin average',
]

ordered = sorted(zip(handles, label), key=lambda hl: desired_order.index(hl[1]))
handles, label = zip(*ordered)


# Customize plot
ax.set_xlabel('CORA centroid', fontsize=12)
ax.set_ylabel('Joint return period (years)', fontsize=12)
ax.set_xticks(cora_ids)
ax.set_xticklabels(labels)
ax.set_ylim(20, 105)  # Adjust as needed
ax.grid(True, alpha=0.3, linestyle='--')
#ax.legend(handles, label, loc='upper center', frameon=True, shadow=True)
ax.legend(handles, label, loc='upper center', 
          bbox_to_anchor=(0.4, 1.0),  # x, y in axes-fraction coords; y>1 pushes it above the plot area
          frameon=True, shadow=True, ncol=1)


#save plot
plt.savefig(f'C:/Users/ro327497/OneDrive - University of Central Florida/Documents/CONUS/NJ/{name}/Wahl approach/{name}_compairon_of_methods_10_yr.png',bbox_inches='tight',pad_inches=0.1)




## Plot correlations for basin average vs correlations in condititonasl 