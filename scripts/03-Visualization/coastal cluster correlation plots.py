#Import relevent packages
import xarray as xr
import geopandas as gpd
from shapely.geometry import MultiLineString, LineString, Point
import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import cartopy.io.img_tiles as cimgt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from matplotlib.lines import Line2D
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.ticker import FormatStrFormatter

#Basin name
name='Mullica'

#AORC and CORA datasets in the watershed of interest (inc. centroids)
aorc_subset = pd.read_csv(f'C:/Users/ro327497/OneDrive - University of Central Florida/Documents/CONUS/NJ/{name}/{name}_aorc_subset.csv')
aorc_centroids_gdf = pd.read_csv(f'C:/Users/ro327497/OneDrive - University of Central Florida/Documents/CONUS/NJ/{name}/{name}_aorc_centroid_lat_lon.csv')

cora_clean = pd.read_csv(f'C:/Users/ro327497/OneDrive - University of Central Florida/Documents/CONUS/NJ/{name}/{name}_cora_clean.csv')
cora_centroids_gdf = pd.read_csv(f'C:/Users/ro327497/OneDrive - University of Central Florida/Documents/CONUS/NJ/{name}/{name}_cora_centroid_lat_lon.csv')


#For rain duration = 20
rain_dur = 24

#Reading in data
dfs = []

for step in range(50, 2400, 50):
    name_file = f'{name}_cor_by_acc_time_{step}'
    df = pd.read_csv(
        f'C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/correlations/{name_file}.csv',
        index_col=False,
        usecols=[0, 1, 2, 5, 15]
    )
    dfs.append(df)

# Concatenate row-wise
data = pd.concat(dfs, axis=0, ignore_index=True)

data = data[data['rain_dur_seq']==24]

#EPY = 5
EPY = 5

# Create a custom colormap
colors = ["#0000ff", "#00ffff", "#ffff00", "#ff7f00", "#ff0000"]  # Define custom colors
custom_cmap = LinearSegmentedColormap.from_list("custom_cmap", colors)  # Create custom colormap

# Adjust color limits
vmin = data.iloc[:,3].min()  # Minimum value in the dataset
vmax = data.iloc[:,3].max()  # Maximum value in the dataset

# Create a ScalarMappable for the colorbar
norm = plt.Normalize(vmin=vmin, vmax=vmax)
sm = plt.cm.ScalarMappable(cmap=custom_cmap, norm=norm)
sm.set_array([])

#Define labels
labels = ['A','B','C','D','E']

#Define terrain
stamen_terrain = cimgt.StadiaMapsTiles(apikey='1dabf9e0-8140-4bfc-9327-2c1ef91e4557', style='stamen_terrain_background', resolution='@2x')


#Create the figure
fig = plt.figure(figsize=(15, 3.5))
 
# Create subplot grid for manual assignment
gs = fig.add_gridspec(1, 5, wspace=0.01)

# Create each subplot individually with the correct projection
axes = []
for i in range(5):  # Create 5 subplots (3×2 grid minus 1)
 row = 0 #i // 2
 col = i #% 2
 ax = fig.add_subplot(gs[row, col], projection=stamen_terrain.crs)
 #ax = fig.add_axes(positions[i], projection=stamen_terrain.crs)
 axes.append(ax)

 cmap = plt.get_cmap('plasma_r')

 # set miami-dade county as map extent
 ax.set_extent([-75, -74, 39, 40], crs=ccrs.Geodetic())

 # Add the Stamen data at zoom level 8.
 ax.add_image(stamen_terrain, 8)
 ax.add_feature(cfeature.COASTLINE)
 ax.add_feature(cfeature.STATES)
 ax.add_feature(cfeature.RIVERS)
 ax.add_feature(cfeature.LAKES)
 ax.add_feature(cfeature.OCEAN)
 ax.add_feature(cfeature.LAND)
    
 # Use scatter to control marker size and transformation
 ax.scatter(aorc_subset.lon, aorc_subset.lat, c=data[data['cora']==i+1].iloc[:, 3], s=1, cmap=custom_cmap, vmin=vmin, vmax=vmax, transform=ccrs.PlateCarree())

 #Add CORA point
 point = cora_centroids_gdf.iloc[i]
 ax.scatter(point.lon, point.lat, color='black', s=20, edgecolors='black', linewidths=1, transform=ccrs.Geodetic())

 mean_val = data[data['cora']==i+1].iloc[:, 4].mean()
 sd_val = data[data['cora']==i+1].iloc[:, 4].std()

 # Add text to the plot
 ax.text(0.95, 0.05, f'# events: {mean_val:.0f}\nsd: {sd_val:.2f}', 
        transform=ax.transAxes,
        fontsize=10,
        verticalalignment='bottom',
        horizontalalignment='right',
        bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
   
 # Add title with label
 ax.set_title(f'CORA centroid {labels[i]}', fontsize=12)
 
# Add colorbar to the figure
cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.7])  # [left, bottom, width, height]
cbar = fig.colorbar(sm, cax=cbar_ax)
cbar.set_label('Correlation', fontsize=12)

# Adjust layout to make room for colorbar
plt.subplots_adjust(right=0.9)

# Display the plot
plt.tight_layout(rect=[0, 0, 0.9, 1])  # Adjust layout but leave space for colorbar

#save plot
plt.savefig(f'C:/Users/ro327497/OneDrive - University of Central Florida/Desktop/conferences/Coastal cluster/anv_{name}_con_{rain_dur}_hr_rain_total_tc_epy_{EPY}_map.png')

plt.show()




##Non-TC
#Different map for each CORA centroid and all AORC points


#Define labels
labels = ['A','B','C','D','E']

#For rain duration = 20
rain_dur = 24

#Reading in data
dfs = []

for step in range(50, 2400, 50):
    name_file = f'{name}_cor_by_acc_time_{step}'
    df = pd.read_csv(
        f'C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/correlations/{name_file}.csv',
        index_col=False,
        usecols=[0, 1, 2, 6, 16]
    )
    dfs.append(df)

# Concatenate row-wise
data = pd.concat(dfs, axis=0, ignore_index=True)

data = data[data['rain_dur_seq']==24]

# Create a custom colormap
colors = ["#0000ff", "#00ffff", "#ffff00", "#ff7f00", "#ff0000"]  # Define custom colors
custom_cmap = LinearSegmentedColormap.from_list("custom_cmap", colors)  # Create custom colormap

# Adjust color limits
vmin = data.iloc[:,3].min()  # Minimum value in the dataset
vmax = data.iloc[:,3].max()  # Maximum value in the dataset

# Create a ScalarMappable for the colorbar
norm = plt.Normalize(vmin=vmin, vmax=vmax)
sm = plt.cm.ScalarMappable(cmap=custom_cmap, norm=norm)
sm.set_array([])

#Define terrain
stamen_terrain = cimgt.StadiaMapsTiles(apikey='1dabf9e0-8140-4bfc-9327-2c1ef91e4557', style='stamen_terrain_background', resolution='@2x')

#Create the figure
fig = plt.figure(figsize=(15, 3.5))
 
# Create subplot grid for manual assignment
gs = fig.add_gridspec(1, 5, wspace=0.01)

# Create each subplot individually with the correct projection
axes = []
for i in range(5):  # Create 5 subplots (3×2 grid minus 1)
 row = 0 #i // 2
 col = i #% 2
 ax = fig.add_subplot(gs[row, col], projection=stamen_terrain.crs)
 #ax = fig.add_axes(positions[i], projection=stamen_terrain.crs)
 axes.append(ax)

 cmap = plt.get_cmap('plasma_r')

 # set miami-dade county as map extent
 ax.set_extent([-75, -74, 39, 40], crs=ccrs.Geodetic())
 # Add the Stamen data at zoom level 8.

 ax.add_image(stamen_terrain, 8)
 ax.add_feature(cfeature.COASTLINE)
 ax.add_feature(cfeature.STATES)
 ax.add_feature(cfeature.RIVERS)
 ax.add_feature(cfeature.LAKES)
 ax.add_feature(cfeature.OCEAN)
 ax.add_feature(cfeature.LAND)
    
 # Use scatter to control marker size and transformation
 ax.scatter(aorc_subset.lon, aorc_subset.lat, c=data[data['cora']==i+1].iloc[:, 3], s=1, cmap=custom_cmap, vmin=vmin, vmax=vmax, transform=ccrs.PlateCarree())

 #Add CORA point
 point = cora_centroids_gdf.iloc[i]
 ax.scatter(point.lon, point.lat, color='black', s=20, transform=ccrs.Geodetic())

 mean_val = data[data['cora']==i+1].iloc[:, 4].mean()
 sd_val = data[data['cora']==i+1].iloc[:, 4].std()

 # Add text to the plot
 ax.text(0.95, 0.05, f'# events: {mean_val:.0f}\nsd: {sd_val:.2f}', 
        transform=ax.transAxes,
        fontsize=10,
        verticalalignment='bottom',
        horizontalalignment='right',
        bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

 # Add title with label
 ax.set_title(f'CORA centroid {labels[i]}', fontsize=12)
 
# Add colorbar to the figure
cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.7])  # [left, bottom, width, height]
cbar = fig.colorbar(sm, cax=cbar_ax)
cbar.set_label('Correlation', fontsize=12)

# Adjust layout to make room for colorbar
plt.subplots_adjust(right=0.9)

# Display the plot
plt.tight_layout(rect=[0, 0, 0.9, 1])  # Adjust layout but leave space for colorbar

#save plot
plt.savefig(f'C:/Users/ro327497/OneDrive - University of Central Florida/Desktop/conferences/Coastal cluster/anv_{name}_con_{rain_dur}_hr_rain_total_non-tc_epy_{EPY}_map.png')

plt.show()


###Conditioned on WL
#Define labels
labels = ['A','B','C','D','E']

##Loop through epy
#Correltations
cols = ['aorc', 'cora', 'cor_1_epy_tc','cor_2_epy_tc','cor_3_epy_tc','cor_4_epy_tc','cor_5_epy_tc', 'n_1_epy_tc', 'n_2_epy_tc', 'n_3_epy_tc', 'n_4_epy_tc', 'n_5_epy_tc']
# Read all three CSVs and concatenate by columns
dfs = []
files = [500, 1000, 1500, 2000, 2300]
for i in files:
   Name =  f'{name}_anv_spatial_acc_time_con_wl_res_detrend_{i}'
   df = pd.read_csv(
        f'C:/Users/ro327497/OneDrive - University of Central Florida/Documents/Paper/prav/results/con_wl/data/{Name}.csv',
        index_col=False, 
        usecols=cols
    )
   dfs.append(df)

# Concatenate by columns (axis=1)
data = pd.concat(dfs, axis=0, ignore_index=True)

#Ensure specified column ordeer
data = data[cols]


#Define terrain
stamen_terrain = cimgt.StadiaMapsTiles(apikey='1dabf9e0-8140-4bfc-9327-2c1ef91e4557', style='stamen_terrain_background', resolution='@2x')


#Loop through different number of events per year
for k in range(1,6):
 data_sub = data.iloc[:,[0,1,k+1, k+6]]
 
 #EPY
 EPY = [1, 2, 3, 4, 5][k-1]
 
 # Create a custom colormap
 colors = ["#0000ff", "#00ffff", "#ffff00", "#ff7f00", "#ff0000"]  # Define custom colors
 custom_cmap = LinearSegmentedColormap.from_list("custom_cmap", colors)  # Create custom colormap

 # Adjust color limits
 vmin = data_sub.iloc[:,2].min()  # Minimum value in the dataset
 vmax = data_sub.iloc[:,2].max()  # Maximum value in the dataset

 # Create a ScalarMappable for the colorbar
 norm = plt.Normalize(vmin=vmin, vmax=vmax)
 sm = plt.cm.ScalarMappable(cmap=custom_cmap, norm=norm)
 sm.set_array([])


 #Create the figure
 fig = plt.figure(figsize=(15, 3.5))
 
 # Create subplot grid for manual assignment
 gs = fig.add_gridspec(1, 5, wspace=0.01)

 # Create each subplot individually with the correct projection
 axes = []
 for i in range(5):  # Create 5 subplots (3×2 grid minus 1)
  row = 0 #i // 2
  col = i #% 2
  ax = fig.add_subplot(gs[row, col], projection=stamen_terrain.crs)
  #ax = fig.add_axes(positions[i], projection=stamen_terrain.crs)
  axes.append(ax)

  cmap = plt.get_cmap('plasma_r')

  # set miami-dade county as map extent
  ax.set_extent([-75, -74, 39, 40], crs=ccrs.Geodetic())
  # Add the Stamen data at zoom level 8.

  ax.add_image(stamen_terrain, 8)
  ax.add_feature(cfeature.COASTLINE)
  ax.add_feature(cfeature.STATES)
  ax.add_feature(cfeature.RIVERS)
  ax.add_feature(cfeature.LAKES)
  ax.add_feature(cfeature.OCEAN)
  ax.add_feature(cfeature.LAND)
    
  # Use scatter to control marker size and transformation
  ax.scatter(aorc_subset.lon, aorc_subset.lat, c=data_sub[data_sub['cora']==i+1].iloc[:, 2], s=1, cmap=custom_cmap, vmin=vmin, vmax=vmax, transform=ccrs.PlateCarree())

  #Add CORA point
  point = cora_centroids_gdf.iloc[[i]]
  ax.scatter(point.lon, point.lat, color='black', s=20, edgecolors='black', linewidths=1, transform=ccrs.Geodetic())

  mean_val = data_sub[data_sub['cora']==i+1].iloc[:, 3].mean()
  sd_val = data_sub[data_sub['cora']==i+1].iloc[:, 3].std()

  # Add text to the plot
  ax.text(0.95, 0.05, f'# events: {mean_val:.0f}\nsd: {sd_val:.2f}', 
          transform=ax.transAxes,
          fontsize=10,
          verticalalignment='bottom',
          horizontalalignment='right',
          bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
  # Add title with label
  ax.set_title(f'CORA centroid {labels[i]}', fontsize=12)
 
 # Add colorbar to the figure
 cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.7])  # [left, bottom, width, height]
 cbar = fig.colorbar(sm, cax=cbar_ax)
 cbar.set_label('Correlation', fontsize=12)

 # Adjust layout to make room for colorbar
 plt.subplots_adjust(right=0.9)

 # Display the plot
 plt.tight_layout(rect=[0, 0, 0.9, 1])  # Adjust layout but leave space for colorbar

 #save plot
 plt.savefig(f'C:/Users/ro327497/OneDrive - University of Central Florida/Desktop/conferences/Coastal cluster/{name}_anv_spatial_24_hr_acc_time_wl_con_wl_tc_epy_{EPY}_map.png')

 plt.show()




##Non-TC
#Different map for each CORA centroid and all AORC points

#Define labels
labels = ['A','B','C','D','E']

##Loop through epy
#Correltations
cols = ['aorc', 'cora', 'cor_1_epy_non_tc','cor_2_epy_non_tc','cor_3_epy_non_tc','cor_4_epy_non_tc','cor_5_epy_non_tc', 'n_1_epy_non_tc', 'n_2_epy_non_tc', 'n_3_epy_non_tc', 'n_4_epy_non_tc', 'n_5_epy_non_tc']
# Read all three CSVs and concatenate by columns
dfs = []
files = [500, 1000, 1500, 2000, 2300]
for i in files:
   Name =  f'{name}_anv_spatial_acc_time_con_wl_res_detrend_{i}'
   df = pd.read_csv(
        f'C:/Users/ro327497/OneDrive - University of Central Florida/Documents/Paper/prav/results/con_wl/data/{Name}.csv',
        index_col=False, 
        usecols=cols
    )
   dfs.append(df)

# Concatenate by columns (axis=1)
data = pd.concat(dfs, axis=0, ignore_index=True)

#Ensure specified column ordeer
data = data[cols]


#Define terrain
stamen_terrain = cimgt.StadiaMapsTiles(apikey='1dabf9e0-8140-4bfc-9327-2c1ef91e4557', style='stamen_terrain_background', resolution='@2x')


#Loop through different number of events per year
for k in range(1,6):
 data_sub = data.iloc[:,[0,1,k+1, k+6]]
 
 #EPY
 EPY = [1, 2, 3, 4, 5][k-1]
 
 # Create a custom colormap
 colors = ["#0000ff", "#00ffff", "#ffff00", "#ff7f00", "#ff0000"]  # Define custom colors
 custom_cmap = LinearSegmentedColormap.from_list("custom_cmap", colors)  # Create custom colormap

 # Adjust color limits
 vmin = data_sub.iloc[:,2].min()  # Minimum value in the dataset
 vmax = data_sub.iloc[:,2].max()  # Maximum value in the dataset

 # Create a ScalarMappable for the colorbar
 norm = plt.Normalize(vmin=vmin, vmax=vmax)
 sm = plt.cm.ScalarMappable(cmap=custom_cmap, norm=norm)
 sm.set_array([])


 #Create the figure
 fig = plt.figure(figsize=(15, 3.5))
 
 # Create subplot grid for manual assignment
 gs = fig.add_gridspec(1, 5, wspace=0.01)

 # Create each subplot individually with the correct projection
 axes = []
 for i in range(5):  # Create 5 subplots (3×2 grid minus 1)
  row = 0 #i // 2
  col = i #% 2
  ax = fig.add_subplot(gs[row, col], projection=stamen_terrain.crs)
  #ax = fig.add_axes(positions[i], projection=stamen_terrain.crs)
  axes.append(ax)

  cmap = plt.get_cmap('plasma_r')

  # set miami-dade county as map extent
  ax.set_extent([-75, -74, 39, 40], crs=ccrs.Geodetic())
  # Add the Stamen data at zoom level 8.

  ax.add_image(stamen_terrain, 8)
  ax.add_feature(cfeature.COASTLINE)
  ax.add_feature(cfeature.STATES)
  ax.add_feature(cfeature.RIVERS)
  ax.add_feature(cfeature.LAKES)
  ax.add_feature(cfeature.OCEAN)
  ax.add_feature(cfeature.LAND)
    
  # Use scatter to control marker size and transformation
  ax.scatter(aorc_subset.lon, aorc_subset.lat, c=data_sub[data_sub['cora']==i+1].iloc[:, 2], s=1, cmap=custom_cmap, vmin=vmin, vmax=vmax, transform=ccrs.PlateCarree())

  #Add CORA point
  point = cora_centroids_gdf.iloc[[i]]
  ax.scatter(point.lon, point.lat, color='black', s=20, edgecolors='black', linewidths=1, transform=ccrs.Geodetic())

  mean_val = data_sub[data_sub['cora']==i+1].iloc[:, 3].mean()
  sd_val = data_sub[data_sub['cora']==i+1].iloc[:, 3].std()

  # Add text to the plot
  ax.text(0.95, 0.05, f'# events: {mean_val:.0f}\nsd: {sd_val:.2f}', 
          transform=ax.transAxes,
          fontsize=10,
          verticalalignment='bottom',
          horizontalalignment='right',
          bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
  # Add title with label
  ax.set_title(f'CORA centroid {labels[i]}', fontsize=12)
 
 # Add colorbar to the figure
 cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.7])  # [left, bottom, width, height]
 cbar = fig.colorbar(sm, cax=cbar_ax)
 cbar.set_label('Correlation', fontsize=12)

 # Adjust layout to make room for colorbar
 plt.subplots_adjust(right=0.9)

 # Display the plot
 plt.tight_layout(rect=[0, 0, 0.9, 1])  # Adjust layout but leave space for colorbar

 #save plot
 plt.savefig(f'C:/Users/ro327497/OneDrive - University of Central Florida/Desktop/conferences/Coastal cluster/{name}_anv_spatial_24_hr_acc_time_wl_con_wl_non_tc_epy_{EPY}_map.png')

 plt.show()

#####NTR not WL
#Import relevent packages
import xarray as xr
import geopandas as gpd
from shapely.geometry import MultiLineString, LineString, Point
import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import cartopy.io.img_tiles as cimgt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from matplotlib.lines import Line2D
from matplotlib.colors import LinearSegmentedColormap

#Basin name
name='Mullica'

#AORC and CORA datasets in the watershed of interest (inc. centroids)
aorc_subset = pd.read_csv(f'C:/Users/ro327497/OneDrive - University of Central Florida/Documents/CONUS/NJ/{name}/{name}_aorc_subset.csv')
aorc_centroids_gdf = pd.read_csv(f'C:/Users/ro327497/OneDrive - University of Central Florida/Documents/CONUS/NJ/{name}/{name}_aorc_centroid_lat_lon.csv')

cora_clean = pd.read_csv(f'C:/Users/ro327497/OneDrive - University of Central Florida/Documents/CONUS/NJ/{name}/{name}_cora_clean.csv')
cora_centroids_gdf = pd.read_csv(f'C:/Users/ro327497/OneDrive - University of Central Florida/Documents/CONUS/NJ/{name}/{name}_cora_centroid_lat_lon.csv')

from shapely import wkt


# Convert WKT strings to Shapely geometries
cora_clean['geometry'] = cora_clean['geometry'].apply(wkt.loads)


# If cora_clean is still a regular DataFrame
cora_clean = gpd.GeoDataFrame(cora_clean, geometry='geometry')
#aorc_subset = gpd.GeoDataFrame(aorc_subset, geometry='geometry')

#Download map
#Different map for each CORA centroid and all AORC points

#Define labels
labels = ['A','B','C','D','E']

##Loop through epy
#Correltations
cols = ['aorc', 'cora', 'cor_est_1_epy_tc','cor_est_2_epy_tc','cor_est_3_epy_tc','cor_est_4_epy_tc','cor_est_5_epy_tc', 'n_1_epy_tc', 'n_2_epy_tc', 'n_3_epy_tc', 'n_4_epy_tc', 'n_5_epy_tc']
# Read all three CSVs and concatenate by columns
dfs = []
files = [500, 1000, 1500, 2000, 2500]
for i in files:
   Name =  f'{name}_anv_spatial_24_acc_time_results_detrend_con_r_with_ntr_{i}'
   df = pd.read_csv(
        f'C:/Users/ro327497/OneDrive - University of Central Florida/Documents/Paper/prav/results/con_ntr/data/{Name}.csv',
        index_col=False, 
        usecols=cols
    )
   dfs.append(df)

# Concatenate by columns (axis=1)
data = pd.concat(dfs, axis=0, ignore_index=True)

#Ensure specified column ordeer
data = data[cols]


#Define terrain
stamen_terrain = cimgt.StadiaMapsTiles(apikey='1dabf9e0-8140-4bfc-9327-2c1ef91e4557', style='stamen_terrain_background', resolution='@2x')


#Loop through different number of events per year
for k in range(1,6):
 data_sub = data.iloc[:,[0,1,k+1, k+6]]
 
 #EPY
 EPY = [1, 2, 3, 4, 5][k-1]
 
 # Create a custom colormap
 colors = ["#0000ff", "#00ffff", "#ffff00", "#ff7f00", "#ff0000"]  # Define custom colors
 custom_cmap = LinearSegmentedColormap.from_list("custom_cmap", colors)  # Create custom colormap

 # Adjust color limits
 vmin = data_sub.iloc[:,2].min()  # Minimum value in the dataset
 vmax = data_sub.iloc[:,2].max()  # Maximum value in the dataset

 # Create a ScalarMappable for the colorbar
 norm = plt.Normalize(vmin=vmin, vmax=vmax)
 sm = plt.cm.ScalarMappable(cmap=custom_cmap, norm=norm)
 sm.set_array([])


 #Create the figure
 fig = plt.figure(figsize=(15, 3.5))
 
 # Create subplot grid for manual assignment
 gs = fig.add_gridspec(1, 5, wspace=0.01)

 # Create each subplot individually with the correct projection
 axes = []
 for i in range(5):  # Create 5 subplots (3×2 grid minus 1)
  row = 0 #i // 2
  col = i #% 2
  ax = fig.add_subplot(gs[row, col], projection=stamen_terrain.crs)
  #ax = fig.add_axes(positions[i], projection=stamen_terrain.crs)
  axes.append(ax)

  cmap = plt.get_cmap('plasma_r')

  # set miami-dade county as map extent
  ax.set_extent([-75, -74, 39, 40], crs=ccrs.Geodetic())
  # Add the Stamen data at zoom level 8.

  ax.add_image(stamen_terrain, 8)
  ax.add_feature(cfeature.COASTLINE)
  ax.add_feature(cfeature.STATES)
  ax.add_feature(cfeature.RIVERS)
  ax.add_feature(cfeature.LAKES)
  ax.add_feature(cfeature.OCEAN)
  ax.add_feature(cfeature.LAND)
    
  # Use scatter to control marker size and transformation
  ax.scatter(aorc_subset.lon, aorc_subset.lat, c=data_sub[data_sub['cora']==i+1].iloc[:, 2], s=1, cmap=custom_cmap, vmin=vmin, vmax=vmax, transform=ccrs.PlateCarree())

  #Add CORA point
  point = cora_centroids_gdf.iloc[[i]]
  ax.scatter(point.lon, point.lat, color='black', s=20, edgecolors='black', linewidths=1, transform=ccrs.Geodetic())

  mean_val = data_sub[data_sub['cora']==i+1].iloc[:, 3].mean()
  sd_val = data_sub[data_sub['cora']==i+1].iloc[:, 3].std()

  # Add text to the plot
  ax.text(0.95, 0.05, f'# events: {mean_val:.0f}\nsd: {sd_val:.2f}', 
          transform=ax.transAxes,
          fontsize=10,
          verticalalignment='bottom',
          horizontalalignment='right',
          bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
  # Add title with label
  ax.set_title(f'CORA centroid {labels[i]}', fontsize=12)
 
 # Add colorbar to the figure
 cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.7])  # [left, bottom, width, height]
 cbar = fig.colorbar(sm, cax=cbar_ax)
 cbar.set_label('Correlation', fontsize=12)

 # Adjust layout to make room for colorbar
 plt.subplots_adjust(right=0.9)

 # Display the plot
 plt.tight_layout(rect=[0, 0, 0.9, 1])  # Adjust layout but leave space for colorbar

 #save plot
 plt.savefig(f'C:/Users/ro327497/OneDrive - University of Central Florida/Desktop/conferences/Coastal cluster/{name}_anv_spatial_24_hr_acc_time_ntr_con_r_tc_epy_{EPY}_map.png')

 plt.show()



# Read all three CSVs and concatenate by columns
cols = ['aorc', 'cora', 'cor_est_1_epy_non_tc','cor_est_2_epy_non_tc','cor_est_3_epy_non_tc','cor_est_4_epy_non_tc','cor_est_5_epy_non_tc', 'n_1_epy_non_tc', 'n_2_epy_non_tc', 'n_3_epy_non_tc', 'n_4_epy_non_tc', 'n_5_epy_non_tc']
dfs = []
files = [500, 1000, 1500, 2000, 2500]
for i in files:
   Name =  f'{name}_anv_spatial_24_acc_time_results_detrend_con_r_with_ntr_{i}'
   df = pd.read_csv(
        f'C:/Users/ro327497/OneDrive - University of Central Florida/Documents/Paper/prav/results/con_ntr/data/{Name}.csv',
        index_col=False, 
        #usecols=[0, 1, 2, 3, 4, 5, 18, 19, 20, 21, 31]
    )
   dfs.append(df)

# Concatenate by columns (axis=1)
data = pd.concat(dfs, axis=0, ignore_index=True)

#Ensuring columns are in specified order
data = data[cols]

#Define terrain
stamen_terrain = cimgt.StadiaMapsTiles(apikey='1dabf9e0-8140-4bfc-9327-2c1ef91e4557', style='stamen_terrain_background', resolution='@2x')


#Loop through different number of events per year
for k in range(1,6):
 data_sub = data.iloc[:,[0,1,k+1, k+6]]
 
 #EPY
 EPY = [1, 2, 3, 4, 5][k-1]
 
 # Create a custom colormap
 colors = ["#0000ff", "#00ffff", "#ffff00", "#ff7f00", "#ff0000"]  # Define custom colors
 custom_cmap = LinearSegmentedColormap.from_list("custom_cmap", colors)  # Create custom colormap

 # Adjust color limits
 vmin = data_sub.iloc[:,2].min()  # Minimum value in the dataset
 vmax = data_sub.iloc[:,2].max()  # Maximum value in the dataset

 # Create a ScalarMappable for the colorbar
 norm = plt.Normalize(vmin=vmin, vmax=vmax)
 sm = plt.cm.ScalarMappable(cmap=custom_cmap, norm=norm)
 sm.set_array([])


 #Create the figure
 fig = plt.figure(figsize=(15, 3.5))
 
 # Create subplot grid for manual assignment
 gs = fig.add_gridspec(1, 5, wspace=0.01)

 # Create each subplot individually with the correct projection
 axes = []
 for i in range(5):  # Create 5 subplots (3×2 grid minus 1)
  row = 0 #i // 2
  col = i #% 2
  ax = fig.add_subplot(gs[row, col], projection=stamen_terrain.crs)
  #ax = fig.add_axes(positions[i], projection=stamen_terrain.crs)
  axes.append(ax)

  cmap = plt.get_cmap('plasma_r')

  # set miami-dade county as map extent
  ax.set_extent([-75, -74, 39, 40], crs=ccrs.Geodetic())
  # Add the Stamen data at zoom level 8.

  ax.add_image(stamen_terrain, 8)
  ax.add_feature(cfeature.COASTLINE)
  ax.add_feature(cfeature.STATES)
  ax.add_feature(cfeature.RIVERS)
  ax.add_feature(cfeature.LAKES)
  ax.add_feature(cfeature.OCEAN)
  ax.add_feature(cfeature.LAND)
    
  # Use scatter to control marker size and transformation
  ax.scatter(aorc_subset.lon, aorc_subset.lat, c=data_sub[data_sub['cora']==i+1].iloc[:, 2], s=1, cmap=custom_cmap, vmin=vmin, vmax=vmax, transform=ccrs.PlateCarree())

  #Add CORA point
  point = cora_centroids_gdf.iloc[[i]]
  ax.scatter(point.lon, point.lat, color='black', s=20, edgecolors='black', linewidths=1, transform=ccrs.Geodetic())

  mean_val = data_sub[data_sub['cora']==i+1].iloc[:, 3].mean()
  sd_val = data_sub[data_sub['cora']==i+1].iloc[:, 3].std()

  # Add text to the plot
  ax.text(0.95, 0.05, f'# events: {mean_val:.0f}\nsd: {sd_val:.2f}', 
          transform=ax.transAxes,
          fontsize=10,
          verticalalignment='bottom',
          horizontalalignment='right',
          bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
  # Add title with label
  ax.set_title(f'CORA centroid {labels[i]}', fontsize=12)
 
 # Add colorbar to the figure
 cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.7])  # [left, bottom, width, height]
 cbar = fig.colorbar(sm, cax=cbar_ax)
 cbar.set_label('Correlation', fontsize=12)

 # Adjust layout to make room for colorbar
 plt.subplots_adjust(right=0.9)

 # Display the plot
 plt.tight_layout(rect=[0, 0, 0.9, 1])  # Adjust layout but leave space for colorbar

 #save plot
 plt.savefig(f'C:/Users/ro327497/OneDrive - University of Central Florida/Desktop/conferences/Coastal cluster/{name}_anv_spatial_24_hr_acc_time_ntr_con_r_non_tc_epy_{EPY}_map.png')

 plt.show()

# Read all three CSVs and concatenate by columns
cols = ['aorc', 'cora', 'cor_1_epy_tc','cor_2_epy_tc','cor_3_epy_tc','cor_4_epy_tc', 'cor_5_epy_tc', 'n_1_epy_tc', 'n_2_epy_tc', 'n_3_epy_tc', 'n_4_epy_tc', 'n_5_epy_tc']
dfs = []
files = [500, 1000, 1500, 2000, 2500]
for i in files:
   Name =  f'{name}_anv_spatial_acc_time_con_ntr_res_detrend_{i}'
   df = pd.read_csv(
        f'C:/Users/ro327497/OneDrive - University of Central Florida/Documents/Paper/prav/results/con_ntr/data/{Name}.csv',
        index_col=False, 
        usecols=cols
    )
   dfs.append(df)

# Concatenate by columns (axis=1)
data = pd.concat(dfs, axis=0, ignore_index=True)

#Ensuring columns are in specified order
data = data[cols]

#Define terrain
stamen_terrain = cimgt.StadiaMapsTiles(apikey='1dabf9e0-8140-4bfc-9327-2c1ef91e4557', style='stamen_terrain_background', resolution='@2x')


#Loop through different number of events per year
for k in range(1,6):
 data_sub = data.iloc[:,[0,1,k+1, k+6]]
 
 #EPY
 EPY = [1, 2, 3, 4, 5][k-1]
 
 # Create a custom colormap
 colors = ["#0000ff", "#00ffff", "#ffff00", "#ff7f00", "#ff0000"]  # Define custom colors
 custom_cmap = LinearSegmentedColormap.from_list("custom_cmap", colors)  # Create custom colormap

 # Adjust color limits
 vmin = data_sub.iloc[:,2].min()  # Minimum value in the dataset
 vmax = data_sub.iloc[:,2].max()  # Maximum value in the dataset

 # Create a ScalarMappable for the colorbar
 norm = plt.Normalize(vmin=vmin, vmax=vmax)
 sm = plt.cm.ScalarMappable(cmap=custom_cmap, norm=norm)
 sm.set_array([])


 #Create the figure
 fig = plt.figure(figsize=(15, 3.5))
 
 # Create subplot grid for manual assignment
 gs = fig.add_gridspec(1, 5, wspace=0.01)

 # Create each subplot individually with the correct projection
 axes = []
 for i in range(5):  # Create 5 subplots (3×2 grid minus 1)
  row = 0 #i // 2
  col = i #% 2
  ax = fig.add_subplot(gs[row, col], projection=stamen_terrain.crs)
  #ax = fig.add_axes(positions[i], projection=stamen_terrain.crs)
  axes.append(ax)

  cmap = plt.get_cmap('plasma_r')

  # set miami-dade county as map extent
  ax.set_extent([-75, -74, 39, 40], crs=ccrs.Geodetic())
  # Add the Stamen data at zoom level 8.

  ax.add_image(stamen_terrain, 8)
  ax.add_feature(cfeature.COASTLINE)
  ax.add_feature(cfeature.STATES)
  ax.add_feature(cfeature.RIVERS)
  ax.add_feature(cfeature.LAKES)
  ax.add_feature(cfeature.OCEAN)
  ax.add_feature(cfeature.LAND)
    
  # Use scatter to control marker size and transformation
  ax.scatter(aorc_subset.lon, aorc_subset.lat, c=data_sub[data_sub['cora']==i+1].iloc[:, 2], s=1, cmap=custom_cmap, vmin=vmin, vmax=vmax, transform=ccrs.PlateCarree())

  #Add CORA point
  point = cora_centroids_gdf.iloc[[i]]
  ax.scatter(point.lon, point.lat, color='black', s=20, edgecolors='black', linewidths=1, transform=ccrs.Geodetic())

  mean_val = data_sub[data_sub['cora']==i+1].iloc[:, 3].mean()
  sd_val = data_sub[data_sub['cora']==i+1].iloc[:, 3].std()

  # Add text to the plot
  ax.text(0.95, 0.05, f'# events: {mean_val:.0f}\nsd: {sd_val:.2f}', 
          transform=ax.transAxes,
          fontsize=10,
          verticalalignment='bottom',
          horizontalalignment='right',
          bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
  # Add title with label
  ax.set_title(f'CORA centroid {labels[i]}', fontsize=12)
 
 # Add colorbar to the figure
 cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.7])  # [left, bottom, width, height]
 cbar = fig.colorbar(sm, cax=cbar_ax)
 cbar.set_label('Correlation', fontsize=12)

 # Adjust layout to make room for colorbar
 plt.subplots_adjust(right=0.9)

 # Display the plot
 plt.tight_layout(rect=[0, 0, 0.9, 1])  # Adjust layout but leave space for colorbar

 #save plot
 plt.savefig(f'C:/Users/ro327497/OneDrive - University of Central Florida/Desktop/conferences/Coastal cluster/{name}_anv_spatial_24_hr_acc_time_ntr_con_ntr_tc_epy_{EPY}_map.png')

 plt.show()
 
 
# Read all three CSVs and concatenate by columns
cols = ['aorc', 'cora', 'cor_1_epy_non_tc','cor_2_epy_non_tc','cor_3_epy_non_tc','cor_4_epy_non_tc', 'cor_5_epy_non_tc', 'n_1_epy_non_tc', 'n_2_epy_non_tc', 'n_3_epy_non_tc', 'n_4_epy_non_tc', 'n_5_epy_non_tc']
dfs = []
files = [500, 1000, 1500, 2000, 2500]
for i in files:
   Name =  f'{name}_anv_spatial_acc_time_con_ntr_res_detrend_{i}'
   df = pd.read_csv(
        f'C:/Users/ro327497/OneDrive - University of Central Florida/Documents/Paper/prav/results/con_ntr/data/{Name}.csv',
        index_col=False, 
        usecols=cols
    )
   dfs.append(df)

# Concatenate by columns (axis=1)
data = pd.concat(dfs, axis=0, ignore_index=True)

#Ensure columns are in specified order
data = data[cols]

#Define terrain
stamen_terrain = cimgt.StadiaMapsTiles(apikey='1dabf9e0-8140-4bfc-9327-2c1ef91e4557', style='stamen_terrain_background', resolution='@2x')


#Loop through different number of events per year
for k in range(1,6):
 data_sub = data.iloc[:,[0,1,k+1, k+6]]
 
 #EPY
 EPY = [1, 2, 3, 4, 5][k-1]
  
 # Create a custom colormap
 colors = ["#0000ff", "#00ffff", "#ffff00", "#ff7f00", "#ff0000"]  # Define custom colors
 custom_cmap = LinearSegmentedColormap.from_list("custom_cmap", colors)  # Create custom colormap

 # Adjust color limits
 vmin = data_sub.iloc[:,2].min()  # Minimum value in the dataset
 vmax = data_sub.iloc[:,2].max()  # Maximum value in the dataset

 # Create a ScalarMappable for the colorbar
 norm = plt.Normalize(vmin=vmin, vmax=vmax)
 sm = plt.cm.ScalarMappable(cmap=custom_cmap, norm=norm)
 sm.set_array([])


 #Create the figure
 fig = plt.figure(figsize=(15, 3.5))
  
 # Create subplot grid for manual assignment
 gs = fig.add_gridspec(1, 5, wspace=0.01)

 # Create each subplot individually with the correct projection
 axes = []
 for i in range(5):  # Create 5 subplots (3×2 grid minus 1)
  row = 0 #i // 2
  col = i #% 2
  ax = fig.add_subplot(gs[row, col], projection=stamen_terrain.crs)
  #ax = fig.add_axes(positions[i], projection=stamen_terrain.crs)
  axes.append(ax)

  cmap = plt.get_cmap('plasma_r')
  # set miami-dade county as map extent
  ax.set_extent([-75, -74, 39, 40], crs=ccrs.Geodetic())
  # Add the Stamen data at zoom level 8.

  ax.add_image(stamen_terrain, 8)
  ax.add_feature(cfeature.COASTLINE)
  ax.add_feature(cfeature.STATES)
  ax.add_feature(cfeature.RIVERS)
  ax.add_feature(cfeature.LAKES)
  ax.add_feature(cfeature.OCEAN)
  ax.add_feature(cfeature.LAND)
    
  # Use scatter to control marker size and transformation
  ax.scatter(aorc_subset.lon, aorc_subset.lat, c=data_sub[data_sub['cora']==i+1].iloc[:, 2], s=1, cmap=custom_cmap, vmin=vmin, vmax=vmax, transform=ccrs.PlateCarree())

  #Add CORA point
  point = cora_centroids_gdf.iloc[[i]]
  ax.scatter(point.lon, point.lat, color='black', s=20, edgecolors='black', linewidths=1, transform=ccrs.Geodetic())

  mean_val = data_sub[data_sub['cora']==i+1].iloc[:, 3].mean()
  sd_val = data_sub[data_sub['cora']==i+1].iloc[:, 3].std()

  # Add text to the plot
  ax.text(0.95, 0.05, f'# events: {mean_val:.0f}\nsd: {sd_val:.2f}', 
          transform=ax.transAxes,
          fontsize=10,
          verticalalignment='bottom',
          horizontalalignment='right',
          bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
  # Add title with label
  ax.set_title(f'CORA centroid {labels[i]}', fontsize=12)
  
 # Add colorbar to the figure
 cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.7])  # [left, bottom, width, height]
 cbar = fig.colorbar(sm, cax=cbar_ax)
 cbar.set_label('Correlation', fontsize=12)

 # Adjust layout to make room for colorbar
 plt.subplots_adjust(right=0.9)
 # Display the plot
 plt.tight_layout(rect=[0, 0, 0.9, 1])  # Adjust layout but leave space for colorbar

 #save plot
 plt.savefig(f'C:/Users/ro327497/OneDrive - University of Central Florida/Desktop/conferences/Coastal cluster/{name}_anv_spatial_24_hr_acc_time_ntr_con_ntr_non_tc_epy_{EPY}_map.png')

plt.show()




























############ Just CORA C for ppt WL vs corresponding NTR ################

# --- WL data tc ---------------------------------------------------------

cols = ['aorc', 'cora', 'cor_r_tc', 'n_r_tc']

dfs = []
steps = list(range(100, 2100, 100)) + [2300]
for step in steps:

    name_file = f'{name}_anv_spatial_24_hr_acc_time_10_yr_r_wl_{step}'

    df = pd.read_csv(
        f'C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/con_wl/data/{name_file}.csv',
        index_col=False,
        usecols=cols
    )

    dfs.append(df)

data_1 = pd.concat(dfs, axis=0, ignore_index=True)

# Keep only 24hr rainfall and CORA=3
data_1 = data_1[data_1['cora'] == 3]


# Keep + rename columns
data_1 = data_1[['aorc', 'cora', 'cor_r_tc', 'n_r_tc']]

data_1.columns = [
    'aorc',
    'cora',
    'cor_5_epy_wl_con_r_tc',
    'n_5_epy_wl_con_r_tc'
]


# --- WL data non-tc -----------------------------------------------------

cols = ['aorc', 'cora', 'cor_r_non_tc', 'n_r_non_tc']

dfs = []
steps = list(range(100, 2100, 100)) + [2300]
for step in steps:

    name_file = f'{name}_anv_spatial_24_hr_acc_time_10_yr_r_wl_{step}'

    df = pd.read_csv(
        f'C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/con_wl/data/{name_file}.csv',
        index_col=False,
        usecols=cols
    )

    dfs.append(df)

data_3 = pd.concat(dfs, axis=0, ignore_index=True)

# Keep only 24hr rainfall and CORA=3
data_3 = data_3[data_3['cora'] == 3]


# Keep + rename columns
data_3 = data_3[['aorc','cora','cor_r_non_tc', 'n_r_non_tc']]

data_3.columns = [
    'aorc',
    'cora',
    'cor_5_epy_wl_con_r_non_tc',
    'n_5_epy_wl_con_r_non_tc'
]


# --- NTR data tc --------------------------------------------------------

cols = ['aorc', 'cora', 'cor_est_5_epy_tc', 'n_5_epy_tc']

dfs = []

files = [500, 1000, 1500, 2000, 2500]

for i in files:

    Name = f'{name}_anv_spatial_24_acc_time_results_detrend_con_r_with_ntr_{i}'

    df = pd.read_csv(
        f'C:/Users/ro327497/OneDrive - University of Central Florida/Documents/Paper/prav/results/con_ntr/data/{Name}.csv',
        index_col=False,
        usecols=cols
    )

    dfs.append(df)

data_2 = pd.concat(dfs, axis=0, ignore_index=True)

# Keep only CORA=3
data_2 = data_2[data_2['cora'] == 3]

# Keep + rename columns
data_2 = data_2[['aorc','cora','cor_est_5_epy_tc', 'n_5_epy_tc']]

data_2.columns = [
    'aorc',
    'cora',
    'cor_5_epy_ntr_con_r_tc',
    'n_5_epy_ntr_con_r_tc'
]


# --- NTR data non-tc ----------------------------------------------------

cols = ['aorc', 'cora', 'cor_est_5_epy_non_tc', 'n_5_epy_non_tc']

dfs = []

files = [500, 1000, 1500, 2000, 2500]

for i in files:

    Name = f'{name}_anv_spatial_24_acc_time_results_detrend_con_r_with_ntr_{i}'

    df = pd.read_csv(
        f'C:/Users/ro327497/OneDrive - University of Central Florida/Documents/Paper/prav/results/con_ntr/data/{Name}.csv',
        index_col=False,
        usecols=cols
    )

    dfs.append(df)

data_4 = pd.concat(dfs, axis=0, ignore_index=True)

# Keep only CORA=3
data_4 = data_4[data_4['cora'] == 3]

# Keep + rename columns
data_4 = data_4[['aorc','cora','cor_est_5_epy_non_tc', 'n_5_epy_non_tc']]

data_4.columns = [
    'aorc',
    'cora',
    'cor_5_epy_ntr_con_r_non_tc',
    'n_5_epy_ntr_con_r_non_tc'
]


# --- Final merged dataframe --------------------------------------------

data = (
    data_1
    .merge(data_2, on=['aorc', 'cora'])
    .merge(data_3, on=['aorc', 'cora'])
    .merge(data_4, on=['aorc', 'cora'])
)


# Final column order
data = data[
    [
        'aorc',
        'cora',
        'cor_5_epy_wl_con_r_tc',
        'n_5_epy_wl_con_r_tc',
        'cor_5_epy_ntr_con_r_tc',
        'n_5_epy_ntr_con_r_tc',
        'cor_5_epy_wl_con_r_non_tc',
        'n_5_epy_wl_con_r_non_tc',
        'cor_5_epy_ntr_con_r_non_tc',
        'n_5_epy_ntr_con_r_non_tc'
    ]
]





#Define terrain
stamen_terrain = cimgt.StadiaMapsTiles(apikey='1dabf9e0-8140-4bfc-9327-2c1ef91e4557', style='stamen_terrain_background', resolution='@2x')

#Row
row_cols = ['cor_5_epy_wl_con_r_tc', 'cor_5_epy_ntr_con_r_tc', 'cor_5_epy_wl_con_r_non_tc', 'cor_5_epy_ntr_con_r_non_tc']
number_cols = ['n_5_epy_wl_con_r_tc', 'n_5_epy_ntr_con_r_tc', 'n_5_epy_wl_con_r_non_tc', 'n_5_epy_ntr_con_r_non_tc']

# Create figure
# Figure
fig = plt.figure(figsize=(8.5, 11))

# Grid
gs = fig.add_gridspec(2, 2, wspace=0.03, hspace=0.02)

# store axes as 2D
axes = [[None for _ in range(2)] for _ in range(4)]

# custom colormap
colors = ["#0000ff", "#00ffff", "#ffff00", "#ff7f00", "#ff0000"]
custom_cmap = LinearSegmentedColormap.from_list("custom_cmap", colors)

# -------------------------
# PLOTTING LOOP
# -------------------------
for i in range(4):

    row = i // 2
    col = i % 2

    ax = fig.add_subplot(gs[row, col], projection=stamen_terrain.crs)
    axes[row][col] = ax

    # subset
    data_sub = data[['aorc', 'cora', row_cols[i], number_cols[i]]]

    # extent
    ax.set_extent([-75, -74, 39, 40], crs=ccrs.PlateCarree())

    # add basemap
    ax.add_image(stamen_terrain, 8)
    ax.add_feature(cfeature.COASTLINE, linewidth=0.5)
    ax.add_feature(cfeature.STATES, linewidth=0.5)
    ax.add_feature(cfeature.RIVERS)
    ax.add_feature(cfeature.LAKES)

    # color limits PER PANEL (you can later change to per-row if needed)
    vmin = data_sub.iloc[:, 2].min()
    vmax = data_sub.iloc[:, 2].max()

    norm = plt.Normalize(vmin=vmin, vmax=vmax)
    sm = plt.cm.ScalarMappable(cmap=custom_cmap, norm=norm)
    sm.set_array([])

    # scatter
    ax.scatter(
        aorc_subset.lon,
        aorc_subset.lat,
        c=data_sub[row_cols[i]],
        s=1,
        cmap=custom_cmap,
        vmin=vmin,
        vmax=vmax,
        transform=ccrs.PlateCarree()
    )

    # CORA point
    point = cora_centroids_gdf.iloc[[2]]
    ax.scatter(
        point.lon,
        point.lat,
        color='black',
        s=20,
        transform=ccrs.PlateCarree()
    )

    # stats box
    mean_val = data_sub.iloc[:, 3].mean()
    sd_val = data_sub.iloc[:, 3].std()

    ax.text(
        0.95, 0.05,
        f'# events: {mean_val:.0f}\nsd: {sd_val:.2f}',
        transform=ax.transAxes,
        fontsize=9,
        ha='right',
        va='bottom',
        bbox=dict(boxstyle='round', facecolor='white', alpha=0.8)
    )

# -------------------------
# ROW-WISE COLORBARS
# -------------------------
row_sm = []

for r in range(2):

    cols = [row_cols[2*r], row_cols[2*r+1]]

    vmin = data[cols].min().min()
    vmax = data[cols].max().max()

    norm = plt.Normalize(vmin, vmax)
    sm = plt.cm.ScalarMappable(cmap=custom_cmap, norm=norm)
    sm.set_array([])

    row_sm.append(sm)

for r in range(2):

    cbar = fig.colorbar(
        row_sm[r],
        ax=axes[r],              # BOTH columns in row
        orientation='vertical',
        fraction=0.035,
        pad=0.02
    )

    cbar.set_label('Correlation', fontsize=12)

# save
plt.savefig(
    f'C:/Users/ro327497/OneDrive - University of Central Florida/Desktop/conferences/Coastal cluster/{name}_ntr_vs_wl_con_r_map.png',
    dpi=300,
    bbox_inches='tight'
)

plt.show()






#Plotting correlation of conditional samples using WT (left column) and NTR (right column)
#Samples conditioned on WLNTR
# --- WL data tc ---------------------------------------------------------

cols = ['aorc', 'cora', 'cor_wl_tc', 'n_wl_tc']

dfs = []
steps = list(range(100, 2100, 100)) + [2300]
for step in steps:

    name_file = f'{name}_anv_spatial_24_hr_acc_time_10_yr_r_wl_{step}'

    df = pd.read_csv(
        f'C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/con_wl/data/{name_file}.csv',
        index_col=False,
        usecols=cols
    )

    dfs.append(df)

data_1 = pd.concat(dfs, axis=0, ignore_index=True)

# Keep only 24hr rainfall and CORA=3
data_1 = data_1[data_1['cora'] == 3]


# Keep + rename columns
data_1 = data_1[['aorc', 'cora', 'cor_wl_tc', 'n_wl_tc']]

data_1.columns = [
    'aorc',
    'cora',
    'cor_5_epy_wl_con_wl_tc',
    'n_5_epy_wl_con_wl_tc'
]


# --- WL data non-tc -----------------------------------------------------

cols = ['aorc', 'cora', 'cor_wl_non_tc', 'n_wl_non_tc']

dfs = []
steps = list(range(100, 2100, 100)) + [2300]
for step in steps:

    name_file = f'{name}_anv_spatial_24_hr_acc_time_10_yr_r_wl_{step}'

    df = pd.read_csv(
        f'C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/con_wl/data/{name_file}.csv',
        index_col=False,
        usecols=cols
    )

    dfs.append(df)

data_3 = pd.concat(dfs, axis=0, ignore_index=True)

# Keep only 24hr rainfall and CORA=3
data_3 = data_3[data_3['cora'] == 3]


# Keep + rename columns
data_3 = data_3[['aorc','cora','cor_wl_non_tc', 'n_wl_non_tc']]

data_3.columns = [
    'aorc',
    'cora',
    'cor_5_epy_wl_con_wl_non_tc',
    'n_5_epy_wl_con_wl_non_tc'
]


# --- NTR data tc --------------------------------------------------------

cols = ['aorc', 'cora', 'cor_5_epy_tc', 'n_5_epy_tc']

dfs = []

files = [500, 1000, 1500, 2000, 2500]

for i in files:

    Name = f'{name}_anv_spatial_acc_time_con_ntr_res_detrend_{i}'

    df = pd.read_csv(
        f'C:/Users/ro327497/OneDrive - University of Central Florida/Documents/Paper/prav/results/con_ntr/data/{Name}.csv',
        index_col=False,
        usecols=cols
    )

    dfs.append(df)

data_2 = pd.concat(dfs, axis=0, ignore_index=True)

# Keep only CORA=3
data_2 = data_2[data_2['cora'] == 3]

# Keep + rename columns
data_2 = data_2[['aorc','cora','cor_5_epy_tc', 'n_5_epy_tc']]

data_2.columns = [
    'aorc',
    'cora',
    'cor_5_epy_ntr_con_ntr_tc',
    'n_5_epy_ntr_con_ntr_tc'
]


# --- NTR data non-tc ----------------------------------------------------

cols = ['aorc', 'cora', 'cor_5_epy_non_tc', 'n_5_epy_non_tc']

dfs = []

files = [500, 1000, 1500, 2000, 2500]

for i in files:

    Name = f'{name}_anv_spatial_acc_time_con_ntr_res_detrend_{i}'

    df = pd.read_csv(
        f'C:/Users/ro327497/OneDrive - University of Central Florida/Documents/Paper/prav/results/con_ntr/data/{Name}.csv',
        index_col=False,
        usecols=cols
    )

    dfs.append(df)

data_4 = pd.concat(dfs, axis=0, ignore_index=True)

# Keep only CORA=3
data_4 = data_4[data_4['cora'] == 3]

# Keep + rename columns
data_4 = data_4[['aorc','cora','cor_5_epy_non_tc', 'n_5_epy_non_tc']]

data_4.columns = [
    'aorc',
    'cora',
    'cor_5_epy_ntr_con_ntr_non_tc',
    'n_5_epy_ntr_con_ntr_non_tc'
]


# --- Final merged dataframe --------------------------------------------

data = (
    data_1
    .merge(data_2, on=['aorc', 'cora'])
    .merge(data_3, on=['aorc', 'cora'])
    .merge(data_4, on=['aorc', 'cora'])
)


# Final column order
data = data[
    [
        'aorc',
        'cora',
        'cor_5_epy_wl_con_wl_tc',
        'n_5_epy_wl_con_wl_tc',
        'cor_5_epy_ntr_con_ntr_tc',
        'n_5_epy_ntr_con_ntr_tc',
        'cor_5_epy_wl_con_wl_non_tc',
        'n_5_epy_wl_con_wl_non_tc',
        'cor_5_epy_ntr_con_ntr_non_tc',
        'n_5_epy_ntr_con_ntr_non_tc'
    ]
]





#Define terrain
stamen_terrain = cimgt.StadiaMapsTiles(apikey='1dabf9e0-8140-4bfc-9327-2c1ef91e4557', style='stamen_terrain_background', resolution='@2x')

#Row
row_cols = ['cor_5_epy_wl_con_wl_tc', 'cor_5_epy_ntr_con_ntr_tc', 'cor_5_epy_wl_con_wl_non_tc', 'cor_5_epy_ntr_con_ntr_non_tc']
number_cols = ['n_5_epy_wl_con_wl_tc', 'n_5_epy_ntr_con_ntr_tc', 'n_5_epy_wl_con_wl_non_tc', 'n_5_epy_ntr_con_ntr_non_tc']

# Create figure
# Figure
fig = plt.figure(figsize=(6,6), constrained_layout=False)

# Grid
gs = fig.add_gridspec(2, 2, wspace=0.03, hspace=0.02)

# store axes as 2D
axes = [[None for _ in range(2)] for _ in range(4)]

# custom colormap
colors = ["#0000ff", "#00ffff", "#ffff00", "#ff7f00", "#ff0000"]
custom_cmap = LinearSegmentedColormap.from_list("custom_cmap", colors)

#Compute min and max of color bar
row_vmin = []
row_vmax = []
for r in range(2):
    cols = [row_cols[2*r], row_cols[2*r+1]]
    row_vmin.append(data[cols].min().min())
    row_vmax.append(data[cols].max().max())
# -------------------------
# PLOTTING LOOP
# -------------------------
for i in range(4):

    row = i // 2
    col = i % 2

    ax = fig.add_subplot(gs[row, col], projection=stamen_terrain.crs)
    axes[row][col] = ax

    # subset
    data_sub = data[['aorc', 'cora', row_cols[i], number_cols[i]]]

    # extent
    ax.set_extent([-75, -74, 39, 40], crs=ccrs.PlateCarree())

    # add basemap
    ax.add_image(stamen_terrain, 8)
    ax.add_feature(cfeature.COASTLINE, linewidth=0.5)
    ax.add_feature(cfeature.STATES, linewidth=0.5)
    ax.add_feature(cfeature.RIVERS)
    ax.add_feature(cfeature.LAKES)

    # color limits PER PANEL (you can later change to per-row if needed)
    vmin = row_vmin[row]
    vmax = row_vmax[row]

    norm = plt.Normalize(vmin=vmin, vmax=vmax)
    sm = plt.cm.ScalarMappable(cmap=custom_cmap, norm=norm)
    sm.set_array([])

    # scatter
    ax.scatter(
        aorc_subset.lon,
        aorc_subset.lat,
        c=data_sub[row_cols[i]],
        s=1,
        cmap=custom_cmap,
        vmin=vmin,
        vmax=vmax,
        transform=ccrs.PlateCarree()
    )

    # CORA point
    point = cora_centroids_gdf.iloc[[2]]
    ax.scatter(
        point.lon,
        point.lat,
        color='black',
        s=20,
        transform=ccrs.PlateCarree()
    )

    # stats box
    mean_val = data_sub.iloc[:, 3].mean()
    sd_val = data_sub.iloc[:, 3].std()

    ax.text(
        0.95, 0.05,
        f'# events: {mean_val:.0f}\nsd: {sd_val:.2f}',
        transform=ax.transAxes,
        fontsize=9,
        ha='right',
        va='bottom',
        bbox=dict(boxstyle='round', facecolor='white', alpha=0.8)
    )


    norm = plt.Normalize(vmin=vmin, vmax=vmax)
    sm = plt.cm.ScalarMappable(cmap=custom_cmap, norm=norm)
    sm.set_array([])

    # scatter
    ax.scatter(
        aorc_subset.lon,
        aorc_subset.lat,
        c=data_sub[row_cols[i]],
        s=1,
        cmap=custom_cmap,
        vmin=vmin,
        vmax=vmax,
        transform=ccrs.PlateCarree()
    )

    # CORA point
    point = cora_centroids_gdf.iloc[[2]]
    ax.scatter(
        point.lon,
        point.lat,
        color='black',
        s=20,
        transform=ccrs.PlateCarree()
    )

    # stats box
    mean_val = data_sub.iloc[:, 3].mean()
    sd_val = data_sub.iloc[:, 3].std()

    ax.text(
        0.95, 0.05,
        f'# events: {mean_val:.0f}\nsd: {sd_val:.2f}',
        transform=ax.transAxes,
        fontsize=9,
        ha='right',
        va='bottom',
        bbox=dict(boxstyle='round', facecolor='white', alpha=0.8)
    )

# -------------------------
# ROW-WISE COLORBARS
# -------------------------

row_titles = ['TC', 'Non-TC']

for r in range(2):
    # get the left axis of each row
    ax = axes[r][0]
    ax.text(
        -0.2, 0.5,          # just outside left edge
        row_titles[r],
        transform=ax.transAxes,
        fontsize=12,
        va='center',
        ha='center',
        rotation=0
    )

row_sm = []

col_titles = ['WL', 'NTR']  # or whatever your columns represent
for c in range(2):
    axes[0][c].set_title(col_titles[c], fontsize=12, pad=3)

# Figure placement
plt.subplots_adjust(
    left=0.13,     # room for row labels
    right=0.88,   # tighter — removes right whitespace
    top=0.98,     
    bottom=0.02   
)

fig.canvas.draw()  # ensures axes positions are finalized before reading them

pos_left = axes[0][0].get_position()
pos_right = axes[0][1].get_position()

x_center = (pos_left.x0 + pos_right.x1) / 2
y_top = pos_left.y1 + 0.05   # tweak offset as needed to clear the WL/NTR titles

fig.text(
    x_center, y_top,
    'Coastal forcing-conditioned',
    fontsize=14,
    ha='center',
    va='bottom'
)


for r in range(2):

    norm = plt.Normalize(row_vmin[r], row_vmax[r])
    sm = plt.cm.ScalarMappable(cmap=custom_cmap, norm=norm)
    sm.set_array([])

    row_sm.append(sm)

for r in range(2):

    cbar = fig.colorbar(
        row_sm[r],
        ax=axes[r],              # BOTH columns in row
        orientation='vertical',
        fraction=0.025,       # was 0.035 — smaller bar
        pad=0.03,            # tighter to panels
        shrink=0.9           # shortens the bar vertically
    )

    cbar.set_label('Correlation', fontsize=12)

# save
plt.savefig(
    f'C:/Users/ro327497/OneDrive - University of Central Florida/Desktop/conferences/Coastal cluster/{name}_ntr_vs_wl_con_wlntr_map.png',
    dpi=300,
    bbox_inches='tight' 
)

plt.show()





































































































#Plotting correlation of conditional samples using WT (left column) and NTR (right column)
#Samples conditioned on rainfall
# --- WL data tc ---------------------------------------------------------

cols = ['aorc', 'cora', 'cor_r_tc', 'n_r_tc']

dfs = []
steps = list(range(100, 2100, 100)) + [2300]
for step in steps:

    name_file = f'{name}_anv_spatial_24_hr_acc_time_10_yr_r_wl_{step}'

    df = pd.read_csv(
        f'C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/con_wl/data/{name_file}.csv',
        index_col=False,
        usecols=cols
    )

    dfs.append(df)

data_1 = pd.concat(dfs, axis=0, ignore_index=True)

# Keep only 24hr rainfall and CORA=3
data_1 = data_1[data_1['cora'] == 3]


# Keep + rename columns
data_1 = data_1[['aorc', 'cora', 'cor_r_tc', 'n_r_tc']]

data_1.columns = [
    'aorc',
    'cora',
    'cor_5_epy_wl_con_r_tc',
    'n_5_epy_wl_con_r_tc'
]


# --- WL data non-tc -----------------------------------------------------

cols = ['aorc', 'cora', 'cor_r_non_tc', 'n_r_non_tc']

dfs = []
steps = list(range(100, 2100, 100)) + [2300]
for step in steps:

    name_file = f'{name}_anv_spatial_24_hr_acc_time_10_yr_r_wl_{step}'

    df = pd.read_csv(
        f'C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/con_wl/data/{name_file}.csv',
        index_col=False,
        usecols=cols
    )

    dfs.append(df)

data_3 = pd.concat(dfs, axis=0, ignore_index=True)

# Keep only 24hr rainfall and CORA=3
data_3 = data_3[data_3['cora'] == 3]


# Keep + rename columns
data_3 = data_3[['aorc','cora','cor_r_non_tc', 'n_r_non_tc']]

data_3.columns = [
    'aorc',
    'cora',
    'cor_5_epy_wl_con_r_non_tc',
    'n_5_epy_wl_con_r_non_tc'
]


# --- NTR data tc --------------------------------------------------------

cols = ['aorc', 'cora', 'cor_est_5_epy_tc', 'n_5_epy_tc']

dfs = []

files = [500, 1000, 1500, 2000, 2500]

for i in files:

    Name = f'{name}_anv_spatial_24_acc_time_results_detrend_con_r_with_ntr_{i}'

    df = pd.read_csv(
        f'C:/Users/ro327497/OneDrive - University of Central Florida/Documents/Paper/prav/results/con_ntr/data/{Name}.csv',
        index_col=False,
        usecols=cols
    )

    dfs.append(df)

data_2 = pd.concat(dfs, axis=0, ignore_index=True)

# Keep only CORA=3
data_2 = data_2[data_2['cora'] == 3]

# Keep + rename columns
data_2 = data_2[['aorc','cora','cor_est_5_epy_tc', 'n_5_epy_tc']]

data_2.columns = [
    'aorc',
    'cora',
    'cor_5_epy_ntr_con_r_tc',
    'n_5_epy_ntr_con_r_tc'
]


# --- NTR data non-tc ----------------------------------------------------

cols = ['aorc', 'cora', 'cor_est_5_epy_non_tc', 'n_5_epy_non_tc']

dfs = []

files = [500, 1000, 1500, 2000, 2500]

for i in files:

    Name = f'{name}_anv_spatial_24_acc_time_results_detrend_con_r_with_ntr_{i}'

    df = pd.read_csv(
        f'C:/Users/ro327497/OneDrive - University of Central Florida/Documents/Paper/prav/results/con_ntr/data/{Name}.csv',
        index_col=False,
        usecols=cols
    )

    dfs.append(df)

data_4 = pd.concat(dfs, axis=0, ignore_index=True)

# Keep only CORA=3
data_4 = data_4[data_4['cora'] == 3]

# Keep + rename columns
data_4 = data_4[['aorc','cora','cor_est_5_epy_non_tc', 'n_5_epy_non_tc']]

data_4.columns = [
    'aorc',
    'cora',
    'cor_5_epy_ntr_con_r_non_tc',
    'n_5_epy_ntr_con_r_non_tc'
]


# --- Final merged dataframe --------------------------------------------

data = (
    data_1
    .merge(data_2, on=['aorc', 'cora'])
    .merge(data_3, on=['aorc', 'cora'])
    .merge(data_4, on=['aorc', 'cora'])
)


# Final column order
data = data[
    [
        'aorc',
        'cora',
        'cor_5_epy_wl_con_r_tc',
        'n_5_epy_wl_con_r_tc',
        'cor_5_epy_ntr_con_r_tc',
        'n_5_epy_ntr_con_r_tc',
        'cor_5_epy_wl_con_r_non_tc',
        'n_5_epy_wl_con_r_non_tc',
        'cor_5_epy_ntr_con_r_non_tc',
        'n_5_epy_ntr_con_r_non_tc'
    ]
]





#Define terrain
stamen_terrain = cimgt.StadiaMapsTiles(apikey='1dabf9e0-8140-4bfc-9327-2c1ef91e4557', style='stamen_terrain_background', resolution='@2x')

#Row
row_cols = ['cor_5_epy_wl_con_r_tc', 'cor_5_epy_ntr_con_r_tc', 'cor_5_epy_wl_con_r_non_tc', 'cor_5_epy_ntr_con_r_non_tc']
number_cols = ['n_5_epy_wl_con_r_tc', 'n_5_epy_ntr_con_r_tc', 'n_5_epy_wl_con_r_non_tc', 'n_5_epy_ntr_con_r_non_tc']

# Create figure
# Figure
fig = plt.figure(figsize=(6,6), constrained_layout=False)

# Grid
gs = fig.add_gridspec(2, 2, wspace=0.03, hspace=0.02)

# store axes as 2D
axes = [[None for _ in range(2)] for _ in range(2)]

# custom colormap
colors = ["#0000ff", "#00ffff", "#ffff00", "#ff7f00", "#ff0000"]
custom_cmap = LinearSegmentedColormap.from_list("custom_cmap", colors)

#Compute min and max of color bar
row_vmin = []
row_vmax = []
for r in range(2):
    cols = [row_cols[2*r], row_cols[2*r+1]]
    row_vmin.append(data[cols].min().min())
    row_vmax.append(data[cols].max().max())
    
# -------------------------
# PLOTTING LOOP
# -------------------------
for i in range(4):

    row = i // 2
    col = i % 2

    ax = fig.add_subplot(gs[row, col], projection=stamen_terrain.crs)
    axes[row][col] = ax

    # subset
    data_sub = data[['aorc', 'cora', row_cols[i], number_cols[i]]]

    # extent
    ax.set_extent([-75, -74, 39, 40], crs=ccrs.PlateCarree())

    # add basemap
    ax.add_image(stamen_terrain, 8)
    ax.add_feature(cfeature.COASTLINE, linewidth=0.5)
    ax.add_feature(cfeature.STATES, linewidth=0.5)
    ax.add_feature(cfeature.RIVERS)
    ax.add_feature(cfeature.LAKES)

    # color limits PER PANEL (you can later change to per-row if needed)
    vmin = row_vmin[row]
    vmax = row_vmax[row]

    norm = plt.Normalize(vmin=vmin, vmax=vmax)
    sm = plt.cm.ScalarMappable(cmap=custom_cmap, norm=norm)
    sm.set_array([])

    # scatter
    ax.scatter(
        aorc_subset.lon,
        aorc_subset.lat,
        c=data_sub[row_cols[i]],
        s=1,
        cmap=custom_cmap,
        vmin=vmin,
        vmax=vmax,
        transform=ccrs.PlateCarree()
    )

    # CORA point
    point = cora_centroids_gdf.iloc[[2]]
    ax.scatter(
        point.lon,
        point.lat,
        color='black',
        s=20,
        transform=ccrs.PlateCarree()
    )

    # stats box
    mean_val = data_sub.iloc[:, 3].mean()
    sd_val = data_sub.iloc[:, 3].std()

    ax.text(
        0.95, 0.05,
        f'# events: {mean_val:.0f}\nsd: {sd_val:.2f}',
        transform=ax.transAxes,
        fontsize=9,
        ha='right',
        va='bottom',
        bbox=dict(boxstyle='round', facecolor='white', alpha=0.8)
    )

# -------------------------
# ROW-WISE COLORBARS
# -------------------------

row_titles = ['TC', 'Non-TC']

for r in range(2):
    # get the left axis of each row
    ax = axes[r][0]
    ax.text(
        -0.2, 0.5,          # just outside left edge
        row_titles[r],
        transform=ax.transAxes,
        fontsize=12,
        va='center',
        ha='center',
        rotation=0
    )

row_sm = []

col_titles = ['WL', 'NTR']  # or whatever your columns represent
for c in range(2):
    axes[0][c].set_title(col_titles[c], fontsize=12,  pad=3)

# Figure placement
plt.subplots_adjust(
    left=0.13,     # room for row labels
    right=0.88,   # tighter — removes right whitespace
    top=0.98,     
    bottom=0.02   
)

fig.canvas.draw()  # ensures axes positions are finalized before reading them

pos_left = axes[0][0].get_position()
pos_right = axes[0][1].get_position()

x_center = (pos_left.x0 + pos_right.x1) / 2
y_top = pos_left.y1 + 0.05   # tweak offset as needed to clear the WL/NTR titles

fig.text(
    x_center, y_top,
    'Rainfall-conditioned',
    fontsize=14,
    ha='center',
    va='bottom'
)

for r in range(2):

    norm = plt.Normalize(row_vmin[r], row_vmax[r])
    sm = plt.cm.ScalarMappable(cmap=custom_cmap, norm=norm)
    sm.set_array([])

    row_sm.append(sm)

for r in range(2):

    cbar = fig.colorbar(
        row_sm[r],
        ax=axes[r],              # BOTH columns in row
        orientation='vertical',
        fraction=0.025,       # was 0.035 — smaller bar
        pad=0.03,            # tighter to panels
        shrink=0.9           # shortens the bar vertically
    )

    cbar.set_label('Correlation', fontsize=12)
    cbar.ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    
# save
plt.savefig(
    f'C:/Users/ro327497/OneDrive - University of Central Florida/Desktop/conferences/Coastal cluster/{name}_ntr_vs_wl_con_r_map.png',
    dpi=300,
    bbox_inches='tight' 
)

plt.show()
