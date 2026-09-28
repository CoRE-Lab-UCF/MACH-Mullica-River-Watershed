#Import relevent packages
import xarray as xr
import geopandas as gpd
from shapely.geometry import MultiLineString, LineString, Point
import pandas as pd
import numpy as np
from shapely import wkt
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

# Convert WKT strings to Shapely geometries
cora_clean['geometry'] = cora_clean['geometry'].apply(wkt.loads)

#Watershed
watersheds = gpd.read_file('C:/Users/ro327497/OneDrive - University of Central Florida/Documents/CONUS/Data/NJ_WMA/depwmas.shp')
watersheds = watersheds.to_crs(epsg=4326)
ws = watersheds[watersheds['WMA_NAME'] == name]



#Define labels
labels = ['A','B','C','D','E']

##Plot exceedence probabilities as proporitons of overall exceedance probability
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


#Define labels
labels = ['A','B','C','D','E']


# Create a custom colormap
colors = ["#0000ff", "#00ffff", "#ffff00", "#ff7f00", "#ff0000"]  # Define custom colors
custom_cmap = LinearSegmentedColormap.from_list("custom_cmap", colors).reversed()  # Create custom colormap

# Adjust color limits
vmin = data['rp_joint'].min()  # Minimum value in the dataset
vmax = data['rp_joint'].max()  # Maximum value in the dataset

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
 ax.scatter(aorc_subset.lon, aorc_subset.lat, c=data[data['cora']==i+1]['rp_joint'], s=1, cmap=custom_cmap, vmin=vmin, vmax=vmax, transform=ccrs.PlateCarree())

 #Add CORA point
 point = cora_centroids_gdf.iloc[[i]]
 ax.scatter(point.lon, point.lat, color='black', s=20, edgecolors='black', linewidths=1, transform=ccrs.Geodetic())

 # Add title with label
 ax.set_title(f'CORA centroid {labels[i]}', fontsize=12)
 
 # Add colorbar to the figure
 cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.7])  # [left, bottom, width, height]
 cbar = fig.colorbar(sm, cax=cbar_ax)
 cbar.set_label('Joint return period', fontsize=12)

 # Adjust layout to make room for colorbar
 plt.subplots_adjust(right=0.9)

 # Display the plot
 plt.tight_layout(rect=[0, 0, 0.9, 1])  # Adjust layout but leave space for colorbar

#save plot
plt.savefig(f'C:/Users/ro327497/OneDrive - University of Central Florida/Desktop/Conferences/Coastal cluster/{name}_10-yr-JRP_map.png')

plt.show()

#How many are greater than 100 a CORA point 3
(data[data['cora']==3]['rp_joint'] > 100).sum()/len(data[data['cora']==3]) #0

# Waht is the maximum return priod for CORA point 3
data[data['cora'] == 3]['rp_joint'].max()

#What percentage are between 30 and 60-years for CORA point 3
((data['cora'] == 3) & (data['rp_joint'] >= 30) & (data['rp_joint'] <= 60)).sum() / len(data[data['cora'] == 3]) #84.5%
 #

#How many are greater than 100 
(data['rp_joint'] > 100).sum() #87
(data['rp_joint'] > 100).sum()/len(data) #0.7%



##Set any above 100 equal to 100 and replot
data['rp_joint'] = data['rp_joint'].clip(upper=100)

#Define labels
labels = ['A','B','C','D','E']


# Create a custom colormap
colors = ["#0000ff", "#00ffff", "#ffff00", "#ff7f00", "#ff0000"]  # Define custom colors
custom_cmap = LinearSegmentedColormap.from_list("custom_cmap", colors).reversed()  # Create custom colormap

# Adjust color limits
vmin = data['rp_joint'].min()  # Minimum value in the dataset
vmax = data['rp_joint'].max()  # Maximum value in the dataset

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
 ax.scatter(aorc_subset.lon, aorc_subset.lat, c=data[data['cora']==i+1]['rp_joint'], s=1, cmap=custom_cmap, vmin=vmin, vmax=vmax, transform=ccrs.PlateCarree())

 #Add CORA point
 point = cora_centroids_gdf.iloc[[i]]
 ax.scatter(point.lon, point.lat, color='black', s=20, edgecolors='black', linewidths=1, transform=ccrs.Geodetic())

 # Add title with label
 ax.set_title(f'CORA centroid {labels[i]}', fontsize=12)
 
 # Add colorbar to the figure
 cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.7])  # [left, bottom, width, height]
 cbar = fig.colorbar(sm, cax=cbar_ax)
 cbar.set_label('Joint return period', fontsize=12)

 # Adjust layout to make room for colorbar
 plt.subplots_adjust(right=0.9)

 # Display the plot
 plt.tight_layout(rect=[0, 0, 0.9, 1])  # Adjust layout but leave space for colorbar

#save plot
plt.savefig(f'C:/Users/ro327497/OneDrive - University of Central Florida/Desktop/Conferences/Coastal cluster/{name}_10-yr-JRP_map_removed.png')

plt.show()



#Finding contributions of sampes to overall annual exceedance probabilities
#Convert AEPs to proportions
cols = [
    'AEP_con_r_tc',
    'AEP_con_r_non_tc',
    'AEP_con_ntr_tc',
    'AEP_con_ntr_non_tc'
]


# Setting correlation value of TC samples assocaited with lower AEP to nan
mask_r_lower = data['AEP_con_r_tc'] <= data['AEP_con_ntr_tc']
data.loc[mask_r_lower, 'cor_r_tc'] = np.nan
data.loc[~mask_r_lower, 'cor_ntr_tc'] = np.nan

# Setting correlation value of non-TC samples assocaited with lower AEP to nan
mask_r_lower = data['AEP_con_r_non_tc'] <= data['AEP_con_ntr_non_tc']
data.loc[mask_r_lower, 'cor_r_non_tc'] = np.nan
data.loc[~mask_r_lower, 'cor_ntr_non_tc'] = np.nan

#Check for where nan's appear
cols = ['cor_r_tc', 'cor_ntr_tc', 'cor_r_non_tc', 'cor_ntr_non_tc']
subset = data[data['cora'] == 3]

result = subset[cols].isna().sum()
print(result)

# Create a custom colormap
colors = ["#0000ff", "#00ffff", "#ffff00", "#ff7f00", "#ff0000"]  # Define custom colors
custom_cmap = LinearSegmentedColormap.from_list("custom_cmap", colors)  # Create custom colormap

# Adjust color limits
vmin = data.iloc[:,2:5].min().min()  # Minimum value in the dataset
vmax = data.iloc[:,2:5].max().max()  # Maximum value in the dataset

# Create a ScalarMappable for the colorbar
norm = plt.Normalize(vmin=vmin, vmax=vmax)
sm = plt.cm.ScalarMappable(cmap=custom_cmap, norm=norm)
sm.set_array([])

#Define terrain
stamen_terrain = cimgt.StadiaMapsTiles(apikey='1dabf9e0-8140-4bfc-9327-2c1ef91e4557', style='stamen_terrain_background', resolution='@2x')


#Loopthrouhg different number of events per year

#Create the figure
fig = plt.figure(figsize=(15, 3.5))
 
# Create subplot grid for manual assignment
gs = fig.add_gridspec(1, 4, wspace=0.01)

#Create labels
labels = ['Rainfall TC','NTR TC', 'Rainfall non-TC', 'NTR non-TC']
cor = cols
# Create each subplot individually with the correct projection
axes = []
for i in range(4):  # Create 5 subplots (3×2 grid minus 1)
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
  ax.scatter(aorc_subset.lon, aorc_subset.lat, c=data[data['cora']==3][cor[i]], s=1, cmap=custom_cmap, vmin=vmin, vmax=vmax, transform=ccrs.PlateCarree())

  #Add CORA point
  point = cora_centroids_gdf.iloc[[2]]
  ax.scatter(point.lon, point.lat, color='black', s=20, edgecolors='black', linewidths=1, transform=ccrs.Geodetic())

  # Add title with label
  ax.set_title(f'{labels[i]}', fontsize=12)
 
# Add colorbar to the figure
cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.7])  # [left, bottom, width, height]
cbar = fig.colorbar(sm, cax=cbar_ax)
cbar.set_label('Correlation', fontsize=12)

# Adjust layout to make room for colorbar
plt.subplots_adjust(right=0.9)

# Display the plot
plt.tight_layout(rect=[0, 0, 0.9, 1])  # Adjust layout but leave space for colorbar

#save plot
plt.savefig(f'C:/Users/ro327497/OneDrive - University of Central Florida/Desktop/Conferences/Coastal cluster/correlations_in_contributing_sample_contributions_{name}.png')

plt.show()

#Correlations are strongest and contribute to JRP
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

# Setting correlation value of TC samples assocaited with lowest correlation to nan
mask_r_lower = data['cor_r_tc'] <= data['cor_ntr_tc']
data.loc[mask_r_lower, 'cor_r_tc'] = np.nan
data.loc[~mask_r_lower, 'cor_ntr_tc'] = np.nan

# Setting correlation value of non-TC samples assocaited with lowest correlation to nan
mask_r_lower = data['cor_r_non_tc'] <= data['cor_ntr_non_tc']
data.loc[mask_r_lower, 'cor_r_non_tc'] = np.nan
data.loc[~mask_r_lower, 'cor_ntr_non_tc'] = np.nan

# Setting correlation value of TC samples assocaited with lower AEP to nan
mask_r_lower = data['AEP_con_r_tc'] <= data['AEP_con_ntr_tc']
data.loc[mask_r_lower, 'cor_r_tc'] = np.nan
data.loc[~mask_r_lower, 'cor_ntr_tc'] = np.nan

# Setting correlation value of non-TC samples assocaited with lower AEP to nan
mask_r_lower = data['AEP_con_r_non_tc'] <= data['AEP_con_ntr_non_tc']
data.loc[mask_r_lower, 'cor_r_non_tc'] = np.nan
data.loc[~mask_r_lower, 'cor_ntr_non_tc'] = np.nan

#Check for where nan's appear
cols = ['cor_r_tc', 'cor_ntr_tc', 'cor_r_non_tc', 'cor_ntr_non_tc']
subset = data[data['cora'] == 3]

result = subset[cols].isna().sum()
print(result)

#Proportion of samples that contribute to AEP and are higher correlation for each sample:
#TC
100 * (len(data[data['cora']==3]['cor_r_tc'].dropna()) + len(data[data['cora']==3]['cor_ntr_tc'].dropna())) / len(data[data['cora']==3]['cora'])
#Non-TC
100 * (len(data[data['cora']==3]['cor_r_non_tc'].dropna()) + len(data[data['cora']==3]['cor_ntr_non_tc'].dropna())) / len(data[data['cora']==3]['cora'])

#Define terrain
stamen_terrain = cimgt.StadiaMapsTiles(apikey='1dabf9e0-8140-4bfc-9327-2c1ef91e4557', style='stamen_terrain_background', resolution='@2x')


#Loopthrouhg different number of events per year

#Create the figure
fig = plt.figure(figsize=(15, 3.5))
 
# Create subplot grid for manual assignment
gs = fig.add_gridspec(1, 4, wspace=0.01)

#Create labels
labels = ['Rainfall TC','NTR TC', 'Rainfall non-TC', 'NTR non-TC']
cor = cols
# Create each subplot individually with the correct projection
axes = []
for i in range(4):  # Create 5 subplots (3×2 grid minus 1)
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
  ax.scatter(aorc_subset.lon, aorc_subset.lat, c=data[data['cora']==3][cor[i]], s=1, cmap=custom_cmap, vmin=vmin, vmax=vmax, transform=ccrs.PlateCarree())

  #Add CORA point
  point = cora_centroids_gdf.iloc[[2]]
  ax.scatter(point.lon, point.lat, color='black', s=20, edgecolors='black', linewidths=1, transform=ccrs.Geodetic())

  # Add title with label
  ax.set_title(f'{labels[i]}', fontsize=12)
 
# Add colorbar to the figure
cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.7])  # [left, bottom, width, height]
cbar = fig.colorbar(sm, cax=cbar_ax)
cbar.set_label('Correlation', fontsize=12)

# Adjust layout to make room for colorbar
plt.subplots_adjust(right=0.9)

# Display the plot
plt.tight_layout(rect=[0, 0, 0.9, 1])  # Adjust layout but leave space for colorbar

#save plot
plt.savefig(f'C:/Users/ro327497/OneDrive - University of Central Florida/Desktop/Conferences/Coastal cluster/correlations_in_contributing_sample_contributions_{name}_strongest.png')

plt.show()


#Correlations are weakest and contribute to JRP
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

# Setting correlation value of TC samples assocaited with stronger correlation to nan
mask_r_lower = data['cor_r_tc'] > data['cor_ntr_tc']
data.loc[mask_r_lower, 'cor_r_tc'] = np.nan
data.loc[~mask_r_lower, 'cor_ntr_tc'] = np.nan

# Setting correlation value of non-TC samples assocaited with stronger correlation to nan
mask_r_lower = data['cor_r_non_tc'] > data['cor_ntr_non_tc']
data.loc[mask_r_lower, 'cor_r_non_tc'] = np.nan
data.loc[~mask_r_lower, 'cor_ntr_non_tc'] = np.nan

# Setting correlation value of TC samples assocaited with lower AEP to nan
mask_r_lower = data['AEP_con_r_tc'] <= data['AEP_con_ntr_tc']
data.loc[mask_r_lower, 'cor_r_tc'] = np.nan
data.loc[~mask_r_lower, 'cor_ntr_tc'] = np.nan

# Setting correlation value of non-TC samples assocaited with lower AEP to nan
mask_r_lower = data['AEP_con_r_non_tc'] <= data['AEP_con_ntr_non_tc']
data.loc[mask_r_lower, 'cor_r_non_tc'] = np.nan
data.loc[~mask_r_lower, 'cor_ntr_non_tc'] = np.nan

#Check for where nan's appear
cols = ['cor_r_tc', 'cor_ntr_tc', 'cor_r_non_tc', 'cor_ntr_non_tc']
subset = data[data['cora'] == 3]

result = subset[cols].isna().sum()
print(result)

#Proportion of samples that contribute to AEP and are weaker correlation for each sample:
#TC (50.3%)
100 * (len(data[data['cora']==3]['cor_r_tc'].dropna()) + len(data[data['cora']==3]['cor_ntr_tc'].dropna())) / len(data[data['cora']==3]['cora'])
#Non-TC (42.46%)
100 * (len(data[data['cora']==3]['cor_r_non_tc'].dropna()) + len(data[data['cora']==3]['cor_ntr_non_tc'].dropna())) / len(data[data['cora']==3]['cora'])


#Check for where nan's appear
cols = ['cor_r_tc', 'cor_ntr_tc', 'cor_r_non_tc', 'cor_ntr_non_tc']
subset = data[data['cora'] == 3]

result = subset[cols].isna().sum()


#Define terrain
stamen_terrain = cimgt.StadiaMapsTiles(apikey='1dabf9e0-8140-4bfc-9327-2c1ef91e4557', style='stamen_terrain_background', resolution='@2x')


#Create the figure
fig = plt.figure(figsize=(15, 3.5))
 
# Create subplot grid for manual assignment
gs = fig.add_gridspec(1, 4, wspace=0.01)

#Create labels
labels = ['Rainfall TC','NTR TC', 'Rainfall non-TC', 'NTR non-TC']
cor = cols
# Create each subplot individually with the correct projection
axes = []
for i in range(4):  # Create 5 subplots (3×2 grid minus 1)
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
  ax.scatter(aorc_subset.lon, aorc_subset.lat, c=subset[cols[i]], s=1, cmap=custom_cmap, vmin=vmin, vmax=vmax, transform=ccrs.PlateCarree())

  #Add CORA point
  point = cora_centroids_gdf.iloc[[2]]
  ax.scatter(point.lon, point.lat, color='black', s=20, edgecolors='black', linewidths=1, transform=ccrs.Geodetic())

  # Add title with label
  ax.set_title(f'{labels[i]}', fontsize=12)
 
# Add colorbar to the figure
cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.7])  # [left, bottom, width, height]
cbar = fig.colorbar(sm, cax=cbar_ax)
cbar.set_label('Correlation', fontsize=12)

# Adjust layout to make room for colorbar
plt.subplots_adjust(right=0.9)

# Display the plot
plt.tight_layout(rect=[0, 0, 0.9, 1])  # Adjust layout but leave space for colorbar

#save plot
plt.savefig(f'C:/Users/ro327497/OneDrive - University of Central Florida/Desktop/Conferences/Coastal cluster/correlations_in_contributing_sample_contributions_{name}_correlation comparisons.png')

plt.show()

#Finding propotion contributions of sampes to overall anual exceedance probabilities

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


# Setting correlation value of TC samples assocaited with lower AEP to nan
mask_r_lower = data['AEP_con_r_tc'] <= data['AEP_con_ntr_tc']
data.loc[mask_r_lower, 'AEP_con_r_tc'] = np.nan
data.loc[~mask_r_lower, 'AEP_con_ntr_tc'] = np.nan

# Setting correlation value of non-TC samples assocaited with lower AEP to nan
mask_r_lower = data['AEP_con_r_non_tc'] <= data['AEP_con_ntr_non_tc']
data.loc[mask_r_lower, 'AEP_con_r_non_tc'] = np.nan
data.loc[~mask_r_lower, 'AEP_con_ntr_non_tc'] = np.nan

#Check for where nan's appear
cols = ['prop_r_tc', 'prop_ntr_tc', 'prop_r_non_tc', 'prop_ntr_non_tc']
subset = data[data['cora'] == 3]

#Proportion of TC and non-TC towards AEP
subset['prop_r_tc'] = subset['AEP_con_r_tc'] / (1-subset['ANEP'])
subset['prop_ntr_tc'] = subset['AEP_con_ntr_tc'] / (1-subset['ANEP'])
subset['prop_r_non_tc']  = subset['AEP_con_r_non_tc'] / (1-subset['ANEP'])
subset['prop_ntr_non_tc']  = subset['AEP_con_ntr_non_tc'] / (1-subset['ANEP'])

# Create a custom colormap
colors = ["#0000ff", "#00ffff", "#ffff00", "#ff7f00", "#ff0000"]  # Define custom colors
custom_cmap = LinearSegmentedColormap.from_list("custom_cmap", colors)  # Create custom colormap

# Adjust color limits
vmin = subset[cols].min().min()  # Minimum value in the dataset
vmax = subset[cols].max().max()  # Maximum value in the dataset

# Create a ScalarMappable for the colorbar
norm = plt.Normalize(vmin=vmin, vmax=vmax)
sm = plt.cm.ScalarMappable(cmap=custom_cmap, norm=norm)
sm.set_array([])

#Define terrain
stamen_terrain = cimgt.StadiaMapsTiles(apikey='1dabf9e0-8140-4bfc-9327-2c1ef91e4557', style='stamen_terrain_background', resolution='@2x')


#Loopthrouhg different number of events per year

#Create the figure
fig = plt.figure(figsize=(15, 3.5))
 
# Create subplot grid for manual assignment
gs = fig.add_gridspec(1, 4, wspace=0.01)

#Create labels
labels = ['Rainfall TC','NTR TC', 'Rainfall non-TC', 'NTR non-TC']
cor = cols
# Create each subplot individually with the correct projection
axes = []
for i in range(4):  # Create 5 subplots (3×2 grid minus 1)
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
  ax.scatter(aorc_subset.lon, aorc_subset.lat, c=subset[cols[i]], s=1, cmap=custom_cmap, vmin=vmin, vmax=vmax, transform=ccrs.PlateCarree())

  #Add CORA point
  point = cora_centroids_gdf.iloc[[2]]
  ax.scatter(point.lon, point.lat, color='black', s=20, edgecolors='black', linewidths=1, transform=ccrs.Geodetic())

  # Add title with label
  ax.set_title(f'{labels[i]}', fontsize=12)
 
# Add colorbar to the figure
cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.7])  # [left, bottom, width, height]
cbar = fig.colorbar(sm, cax=cbar_ax)
cbar.set_label('Proportion of AEP', fontsize=12)

# Adjust layout to make room for colorbar
plt.subplots_adjust(right=0.9)

# Display the plot
plt.tight_layout(rect=[0, 0, 0.9, 1])  # Adjust layout but leave space for colorbar

#save plot
plt.savefig(f'C:/Users/ro327497/OneDrive - University of Central Florida/Desktop/Conferences/Coastal cluster/proportion_in_contributing_sample_contributions_{name}.png')

plt.show()

##relationship between copula family and AEP

#Copulas (TC)
#Dictionary mapping VineCopula package copula numbers to names
vine_copula_dict = {
    0: "Independence",
    1: "Gaussian",
    2: "Student t",
    3: "Clayton",
    4: "Gumbel",
    5: "Frank",
    6: "Joe",
    7: "BB1 (Clayton-Gumbel)",
    8: "BB6 (Joe-Gumbel)",
    9: "BB7 (Joe-Clayton)",
    10: "BB8 (Joe-Frank)",
    13: "Rot. Clayton (180°)",
    14: "Rot. Gumbel (180°)",
    16: "Rot. Joe (180°)",
    17: "Rot. BB1 (180°)",
    18: "Rot. BB6 (180°)",
    19: "Rot. BB7 (180°)",
    20: "Rot. BB8 (180°)",
    23: "Rot. Clayton (90°)",
    24: "Rot. Gumbel (90°)",
    26: "Rot. Joe (90°)",
    27: "Rot. BB1 (90°)",
    28: "Rot. BB6 (90°)",
    29: "Rot. BB7 (90°)",
    30: "Rot. BB8 (90°)",
    33: "Rot. Clayton (270°)",
    34: "Rot. Gumbel (270°)",
    36: "Rot. Joe (270°)",
    37: "Rot. BB1 (270°)",
    38: "Rot. BB6 (270°)",
    39: "Rot. BB7 (270°)",
    40: "Rot. BB8 (270°)",
    104: "Tawn type 1",
    114: "Rot. Tawn type 1 (180°)",
    124: "Rot. Tawn type 1 (90°)",
    134: "Rot. Tawn type 1 (270°)",
    204: "Tawn type 2",
    214: "Rot. Tawn type 2 (180°)",
    224: "Rot. Tawn type 2 (90°)",
    234: "Rot. Tawn type 2 (270°)"
}

vine_copula_font_dict = {
    0: "Independence",
    1: "Gaussian",
    2:  r"$\mathbf{Student t}$",
    3: "Clayton",
    4:  r"$\mathbf{Gumbel}$",
    5:  "Frank",
    6:  r"$\mathbf{Joe}$",
    7:  r"$\mathbf{BB1 (Clayton-Gumbel)}$",
    8:  r"$\mathbf{BB6 (Joe-Gumbel)}$",
    9:  r"$\mathbf{BB7 (Joe-Clayton)}$",
    10: r"$\mathbf{BB8 (Joe-Frank)}$",
    13: r"$\mathbf{Rot. Clayton (180°)}$",
    14: "Rot. Gumbel (180°)",
    16: "Rot. Joe (180°)",
    17: r"$\mathbf{Rot. BB1 (180°)}$",
    18: "Rot. BB6 (180°)",
    19: r"$\mathbf{Rot. BB7 (180°)}$",
    20: "Rot. BB8 (180°)",
    23: "Rot. Clayton (90°)",
    24: "Rot. Gumbel (90°)",
    26: "Rot. Joe (90°)",
    27: "Rot. BB1 (90°)",
    28: "Rot. BB6 (90°)",
    29: "Rot. BB7 (90°)",
    30: "Rot. BB8 (90°)",
    33: "Rot. Clayton (270°)",
    34: "Rot. Gumbel (20°)",
    36: "Rot. Joe (270°)",
    37: "Rot. BB1 (270°)",
    38: "Rot. BB6 (270°)",
    39: "Rot. BB7 (270°)",
    40: "Rot. BB8 (270°)",
    104: r"$\mathbf{Tawn\ type\ 1}$",
    114: "Rot. Tawn type 1 (180°)",
    124: "Rot. Tawn type 1 (90°)",
    134: "Rot. Tawn type 1 (270°)",
    204: r"$\mathbf{Tawn type 2}$",
    214: "Rot. Tawn type 2 (180°)",
    224: "Rot. Tawn type 2 (90°)",
    234: "Rot. Tawn type 2 (270°)"
}

##Copula TC ntr
#Define labels
labels = ['A','B','C','D','E']

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

# Set joint return periods greater than 100 to 100
data['rp_joint'] = data['rp_joint'].clip(upper=100)

#Remove the rows of nans
# Get boolean mask of NaN rows
#nan_mask = data['cop_ntr_tc_family'].isna()
# View the rows with NaN
#data[nan_mask]
# Or just get the indices
#data[nan_mask].index

#data = data.dropna(subset=['cop_wl_tc_family'])

# Setting correlation value of TC samples assocaited with lower AEP to nan
mask_r_lower = data['AEP_con_r_tc'] <= data['AEP_con_ntr_tc']
data.loc[mask_r_lower, 'cop_r_tc_family'] = np.nan
data.loc[~mask_r_lower, 'cop_ntr_tc_family'] = np.nan

# Setting correlation value of non-TC samples assocaited with lower AEP to nan
mask_r_lower = data['AEP_con_r_non_tc'] <= data['AEP_con_ntr_non_tc']
data.loc[mask_r_lower, 'cop_r_non_tc_family'] = np.nan
data.loc[~mask_r_lower, 'cop_ntr_non_tc_family'] = np.nan

#Proportion of sites where 
subset = data[data['cora'] == 3]['cop_ntr_tc_family']
proportion = subset.isin([14]).sum() / len(subset)

#Copulas
copula = data[data['cora']==3]['cop_r_tc_family'].unique()

# Create a custom colormap
colors = ["#0000ff", "#00ffff", "#ffff00", "#ff7f00", "#ff0000"]  # Define custom colors
custom_cmap = LinearSegmentedColormap.from_list("custom_cmap", colors)  # Create custom colormap

# Adjust color limits
vmin = data['rp_joint'].min()  # Minimum value in the dataset
vmax = data['rp_joint'].max()  # Maximum value in the dataset

# Create a ScalarMappable for the colorbar
norm = plt.Normalize(vmin=vmin, vmax=vmax)
sm = plt.cm.ScalarMappable(cmap=custom_cmap, norm=norm)
sm.set_array([])

#Define terrain
stamen_terrain = cimgt.StadiaMapsTiles(apikey='1dabf9e0-8140-4bfc-9327-2c1ef91e4557', style='stamen_terrain_background', resolution='@2x')


#Loopthrouhg different number of events per year

#Create the figure
fig = plt.figure(figsize=(15, 3.5))
 
# Create subplot grid for manual assignment
gs = fig.add_gridspec(1, len(copula), wspace=0.01)



# Create each subplot individually with the correct projection
axes = []
for i in range(len(copula)):  # Create 5 subplots (3×2 grid minus 1)
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
    
  #Set set to just centroid D
  data_subset = data[data['cora']==3]
  data_subset = data_subset.reset_index(drop=True)
  
  #Only selecting sites where copula in sample conditioned on wl is the ith on the list
  mask_cop = data_subset['cop_r_tc_family'] == copula[i]
   
  #Get lebel by mapping copula code to name
  copula_label = vine_copula_dict.get(copula[i], str(copula[i]))
  copula_label_font = vine_copula_font_dict.get(copula[i], str(copula[i]))
  
  # Get the indices where mask_cop is True
  indices_to_plot = data_subset[mask_cop].index
    
   # Use scatter to control marker size and transformation
  ax.scatter(aorc_subset.lon[indices_to_plot], aorc_subset.lat[indices_to_plot], c=data_subset[mask_cop]['rp_joint'], s=1, cmap=custom_cmap, vmin=vmin, vmax=vmax, transform=ccrs.PlateCarree())

  #Add CORA point
  point = cora_centroids_gdf.iloc[[2]]
  ax.scatter(point.lon, point.lat, color='black', s=20, edgecolors='black', linewidths=1, transform=ccrs.Geodetic())

  # Add title with label
  ax.set_title(copula_label_font, fontsize=10)
 
# Add colorbar to the figure
cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.7])  # [left, bottom, width, height]
cbar = fig.colorbar(sm, cax=cbar_ax)
cbar.set_label('Joint return period', fontsize=12)

# Adjust layout to make room for colorbar
plt.subplots_adjust(right=0.9)

# Display the plot
plt.tight_layout(rect=[0, 0, 0.9, 1])  # Adjust layout but leave space for colorbar

#save plot
plt.savefig(f'C:/Users/ro327497/OneDrive - University of Central Florida/Desktop/Conferences/Coastal cluster/{Name}_copula_r_sample_contributions.png')

plt.show()



##Copula TC ntr
#Define labels
labels = ['A','B','C','D','E']

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

# Set joint return periods greater than 100 to 100
data['rp_joint'] = data['rp_joint'].clip(upper=100)

#Remove the rows of nans
# Get boolean mask of NaN rows
#nan_mask = data['cop_ntr_tc_family'].isna()
# View the rows with NaN
#data[nan_mask]
# Or just get the indices
#data[nan_mask].index

#data = data.dropna(subset=['cop_wl_tc_family'])

#Copulas
copula = data[data['cora']==3]['cop_ntr_tc_family'].unique()

# Setting correlation value of TC samples assocaited with lower AEP to nan
mask_r_lower = data['AEP_con_r_tc'] <= data['AEP_con_ntr_tc']
data.loc[mask_r_lower, 'cop_r_tc_family'] = np.nan
data.loc[~mask_r_lower, 'cop_ntr_tc_family'] = np.nan

# Setting correlation value of non-TC samples assocaited with lower AEP to nan
mask_r_lower = data['AEP_con_r_non_tc'] <= data['AEP_con_ntr_non_tc']
data.loc[mask_r_lower, 'cop_r_non_tc_family'] = np.nan
data.loc[~mask_r_lower, 'cop_ntr_non_tc_family'] = np.nan


# Create a custom colormap
colors = ["#0000ff", "#00ffff", "#ffff00", "#ff7f00", "#ff0000"]  # Define custom colors
custom_cmap = LinearSegmentedColormap.from_list("custom_cmap", colors)  # Create custom colormap

# Adjust color limits
vmin = data['rp_joint'].min()  # Minimum value in the dataset
vmax = data['rp_joint'].max()  # Maximum value in the dataset

# Create a ScalarMappable for the colorbar
norm = plt.Normalize(vmin=vmin, vmax=vmax)
sm = plt.cm.ScalarMappable(cmap=custom_cmap, norm=norm)
sm.set_array([])

#Define terrain
stamen_terrain = cimgt.StadiaMapsTiles(apikey='1dabf9e0-8140-4bfc-9327-2c1ef91e4557', style='stamen_terrain_background', resolution='@2x')


#Loopthrouhg different number of events per year

#Create the figure
fig = plt.figure(figsize=(15, 3.5))
 
# Create subplot grid for manual assignment
gs = fig.add_gridspec(1, len(copula), wspace=0.01)



# Create each subplot individually with the correct projection
axes = []
for i in range(len(copula)):  # Create 5 subplots (3×2 grid minus 1)
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
    
  #Add watershed bounrary
  #ws.boundary.plot(ax=ax, edgecolor='black', linewidth=1, label='watershed',  transform=ccrs.Geodetic())
  
  #Set set to just centroid D
  data_subset = data[data['cora']==3]
  data_subset = data_subset.reset_index(drop=True)
  
  #Only selecting sites where copula in sample conditioned on wl is the ith on the list
  mask_cop = data_subset['cop_ntr_tc_family'] == copula[i]
    
  #Get the name of copula for the label
  cop_label = vine_copula_dict.get(copula[i], str(copula[i]))
  copula_label_font = vine_copula_font_dict.get(copula[i], str(copula[i]))
  
  # Get the indices where mask_cop is True
  indices_to_plot = data_subset[mask_cop].index
    
   # Use scatter to control marker size and transformation
  ax.scatter(aorc_subset.lon[indices_to_plot], aorc_subset.lat[indices_to_plot], c=data_subset[mask_cop]['rp_joint'], s=1, cmap=custom_cmap, vmin=vmin, vmax=vmax, transform=ccrs.PlateCarree())

  #Add CORA point
  point = cora_centroids_gdf.iloc[[2]]
  ax.scatter(point.lon, point.lat, color='black', s=20, edgecolors='black', linewidths=1, transform=ccrs.Geodetic())

  # Add title with label
  ax.set_title(copula_label_font, fontsize=10)
 
# Add colorbar to the figure
cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.7])  # [left, bottom, width, height]
cbar = fig.colorbar(sm, cax=cbar_ax)
cbar.set_label('Joint return period', fontsize=12)

# Adjust layout to make room for colorbar
plt.subplots_adjust(right=0.9)

# Display the plot
plt.tight_layout(rect=[0, 0, 0.9, 1])  # Adjust layout but leave space for colorbar

#save plot
plt.savefig(f'C:/Users/ro327497/OneDrive - University of Central Florida/Desktop/Conferences/Coastal cluster/{Name}_copula_ntr_sample_contributions.png')

plt.show()



####All copula families on same plot
##TC
#Define labels
labels = ['A','B','C','D','E']

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

# Set joint return periods greater than 100 to 100
data['rp_joint'] = data['rp_joint'].clip(upper=100)

#Remove the rows of nans
# Get boolean mask of NaN rows
#nan_mask = data['cop_ntr_tc_family'].isna()
# View the rows with NaN
#data[nan_mask]
# Or just get the indices
#data[nan_mask].index

#data = data.dropna(subset=['cop_wl_tc_family'])

# Setting correlation value of TC samples assocaited with lower AEP to nan
mask_r_lower = data['AEP_con_r_tc'] <= data['AEP_con_ntr_tc']
data.loc[mask_r_lower, 'cop_r_tc_family'] = np.nan
data.loc[~mask_r_lower, 'cop_ntr_tc_family'] = np.nan

# Setting correlation value of non-TC samples assocaited with lower AEP to nan
mask_r_lower = data['AEP_con_r_non_tc'] <= data['AEP_con_ntr_non_tc']
data.loc[mask_r_lower, 'cop_r_non_tc_family'] = np.nan
data.loc[~mask_r_lower, 'cop_ntr_non_tc_family'] = np.nan

#Copulas
copula_all_unique = (data[data['cora']==3]['cop_r_tc_family'].dropna().unique(),data[data['cora']==3]['cop_ntr_tc_family'].dropna().unique())
copula_all_unique = np.concatenate(copula_all_unique)

sample = (
    ['cop_r_tc_family'] * len(data.loc[data['cora'] == 3, 'cop_r_tc_family'].dropna().unique()) +
    ['cop_ntr_tc_family'] * len(data.loc[data['cora'] == 3, 'cop_ntr_tc_family'].dropna().unique())
)

# Create a custom colormap
colors = ["#0000ff", "#00ffff", "#ffff00", "#ff7f00", "#ff0000"]  # Define custom colors
custom_cmap = LinearSegmentedColormap.from_list("custom_cmap", colors)  # Create custom colormap

# Adjust color limits
vmin = data['rp_joint'].min()  # Minimum value in the dataset
vmax = data['rp_joint'].max()  # Maximum value in the dataset

# Create a ScalarMappable for the colorbar
norm = plt.Normalize(vmin=vmin, vmax=vmax)
sm = plt.cm.ScalarMappable(cmap=custom_cmap, norm=norm)
sm.set_array([])

#Define terrain
stamen_terrain = cimgt.StadiaMapsTiles(apikey='1dabf9e0-8140-4bfc-9327-2c1ef91e4557', style='stamen_terrain_background', resolution='@2x')


#Loopthrouhg different number of events per year

#Create the figure
fig = plt.figure(figsize=(15, 3.5))
 
# Create subplot grid for manual assignment
gs = fig.add_gridspec(1, len(copula_all_unique), wspace=0.01)



# Create each subplot individually with the correct projection
axes = []
for i in range(len(copula_all_unique)):  # Create 5 subplots (3×2 grid minus 1)
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
    
  #Add watershed bounrary
  #ws.boundary.plot(ax=ax, edgecolor='black', linewidth=1, label='watershed',  transform=ccrs.Geodetic())
  
  #Set set to just centroid D
  data_subset = data[data['cora']==3]
  data_subset = data_subset.reset_index(drop=True)
  
  #Only selecting sites where copula in sample conditioned on wl is the ith on the list
  mask_cop = data_subset[sample[i]] == copula_all_unique[i]
    
  #Get the name of copula for the label
  cop_label = vine_copula_dict.get(copula_all_unique[i], str(copula_all_unique[i]))
  copula_label_font = vine_copula_font_dict.get(copula_all_unique[i], str(copula_all_unique[i]))
  
  # Get the indices where mask_cop is True
  indices_to_plot = data_subset[mask_cop].index
    
   # Use scatter to control marker size and transformation
  ax.scatter(aorc_subset.lon[indices_to_plot], aorc_subset.lat[indices_to_plot], c=data_subset[mask_cop]['rp_joint'], s=1, cmap=custom_cmap, vmin=vmin, vmax=vmax, transform=ccrs.PlateCarree())

  #Add CORA point
  point = cora_centroids_gdf.iloc[[2]]
  ax.scatter(point.lon, point.lat, color='black', s=20, edgecolors='black', linewidths=1, transform=ccrs.Geodetic())

  # Add title with label
  ax.set_title(copula_label_font, fontsize=10)
 
# Add colorbar to the figure
cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.7])  # [left, bottom, width, height]
cbar = fig.colorbar(sm, cax=cbar_ax)
cbar.set_label('Joint return period', fontsize=12)

# Adjust layout to make room for colorbar
plt.subplots_adjust(right=0.9)

# Display the plot
plt.tight_layout(rect=[0, 0, 0.9, 1])  # Adjust layout but leave space for colorbar

#save plot
plt.savefig(f'C:/Users/ro327497/OneDrive - University of Central Florida/Desktop/Conferences/Coastal cluster/{name}_copula_tc_single_plot.png')

plt.show()


##Non-TC
#Define labels
labels = ['A','B','C','D','E']

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

# Set joint return periods greater than 100 to 100
data['rp_joint'] = data['rp_joint'].clip(upper=100)

#Remove the rows of nans
# Get boolean mask of NaN rows
#nan_mask = data['cop_ntr_tc_family'].isna()
# View the rows with NaN
#data[nan_mask]
# Or just get the indices
#data[nan_mask].index

#data = data.dropna(subset=['cop_wl_tc_family'])

# Setting correlation value of TC samples assocaited with lower AEP to nan
mask_r_lower = data['AEP_con_r_tc'] <= data['AEP_con_ntr_tc']
data.loc[mask_r_lower, 'cop_r_tc_family'] = np.nan
data.loc[~mask_r_lower, 'cop_ntr_tc_family'] = np.nan

# Setting correlation value of non-TC samples assocaited with lower AEP to nan
mask_r_lower = data['AEP_con_r_non_tc'] <= data['AEP_con_ntr_non_tc']
data.loc[mask_r_lower, 'cop_r_non_tc_family'] = np.nan
data.loc[~mask_r_lower, 'cop_ntr_non_tc_family'] = np.nan

#Copulas
copula_all_unique = (data[data['cora']==3]['cop_r_non_tc_family'].dropna().unique(),data[data['cora']==3]['cop_ntr_non_tc_family'].dropna().unique())
copula_all_unique = np.concatenate(copula_all_unique)

sample = (
    ['cop_r_non_tc_family'] * len(data.loc[data['cora'] == 3, 'cop_r_non_tc_family'].dropna().unique()) +
    ['cop_ntr_non_tc_family'] * len(data.loc[data['cora'] == 3, 'cop_ntr_non_tc_family'].dropna().unique())
)

# Create a custom colormap
colors = ["#0000ff", "#00ffff", "#ffff00", "#ff7f00", "#ff0000"]  # Define custom colors
custom_cmap = LinearSegmentedColormap.from_list("custom_cmap", colors)  # Create custom colormap

# Adjust color limits
vmin = data['rp_joint'].min()  # Minimum value in the dataset
vmax = data['rp_joint'].max()  # Maximum value in the dataset

# Create a ScalarMappable for the colorbar
norm = plt.Normalize(vmin=vmin, vmax=vmax)
sm = plt.cm.ScalarMappable(cmap=custom_cmap, norm=norm)
sm.set_array([])

#Define terrain
stamen_terrain = cimgt.StadiaMapsTiles(apikey='1dabf9e0-8140-4bfc-9327-2c1ef91e4557', style='stamen_terrain_background', resolution='@2x')


#Loopthrouhg different number of events per year

#Create the figure
fig = plt.figure(figsize=(15, 3.5))
 
# Create subplot grid for manual assignment
gs = fig.add_gridspec(1, len(copula_all_unique), wspace=0.01)



# Create each subplot individually with the correct projection
axes = []
for i in range(len(copula_all_unique)):  # Create 5 subplots (3×2 grid minus 1)
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
    
  #Add watershed bounrary
  #ws.boundary.plot(ax=ax, edgecolor='black', linewidth=1, label='watershed',  transform=ccrs.Geodetic())
  
  #Set set to just centroid D
  data_subset = data[data['cora']==3]
  data_subset = data_subset.reset_index(drop=True)
  
  #Only selecting sites where copula in sample conditioned on wl is the ith on the list
  mask_cop = data_subset[sample[i]] == copula_all_unique[i]
    
  #Get the name of copula for the label
  cop_label = vine_copula_dict.get(copula_all_unique[i], str(copula_all_unique[i]))
  copula_label_font = vine_copula_font_dict.get(copula_all_unique[i], str(copula_all_unique[i]))
  
  # Get the indices where mask_cop is True
  indices_to_plot = data_subset[mask_cop].index
    
   # Use scatter to control marker size and transformation
  ax.scatter(aorc_subset.lon[indices_to_plot], aorc_subset.lat[indices_to_plot], c=data_subset[mask_cop]['rp_joint'], s=1, cmap=custom_cmap, vmin=vmin, vmax=vmax, transform=ccrs.PlateCarree())

  #Add CORA point
  point = cora_centroids_gdf.iloc[[2]]
  ax.scatter(point.lon, point.lat, color='black', s=20, edgecolors='black', linewidths=1, transform=ccrs.Geodetic())

  # Add title with label
  ax.set_title(copula_label_font, fontsize=10)
 
# Add colorbar to the figure
cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.7])  # [left, bottom, width, height]
cbar = fig.colorbar(sm, cax=cbar_ax)
cbar.set_label('Joint return period', fontsize=12)

# Adjust layout to make room for colorbar
plt.subplots_adjust(right=0.9)

# Display the plot
plt.tight_layout(rect=[0, 0, 0.9, 1])  # Adjust layout but leave space for colorbar

#save plot
plt.savefig(f'C:/Users/ro327497/OneDrive - University of Central Florida/Desktop/Conferences/Coastal cluster/{name}_copula_non_tc_single_plot.png')

plt.show()


##Average joint return period of points where either copula has UTD 



#Define labels
labels = ['A','B','C','D','E']

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

# Set joint return periods greater than 100 to 100
data['rp_joint'] = data['rp_joint'].clip(upper=100)

#Remove the rows of nans

# Setting correlation value of TC samples assocaited with lower AEP to nan
mask_r_lower = data['AEP_con_r_tc'] <= data['AEP_con_ntr_tc']
data.loc[mask_r_lower, 'cop_r_tc_family'] = np.nan
data.loc[~mask_r_lower, 'cop_ntr_tc_family'] = np.nan

# Setting correlation value of non-TC samples assocaited with lower AEP to nan
mask_r_lower = data['AEP_con_r_non_tc'] <= data['AEP_con_ntr_non_tc']
data.loc[mask_r_lower, 'cop_r_non_tc_family'] = np.nan
data.loc[~mask_r_lower, 'cop_ntr_non_tc_family'] = np.nan


# Copulas with upper tail dependence
cop_utd = [2, 4, 6, 7, 8, 9, 10, 13, 17, 19, 104, 204]

mask_utd = (
    data['cop_r_tc_family'].isin(cop_utd)
    | data['cop_r_non_tc_family'].isin(cop_utd)
    | data['cop_ntr_tc_family'].isin(cop_utd)
    | data['cop_ntr_non_tc_family'].isin(cop_utd)
)

print(data.loc[mask_utd, 'rp_joint'].mean())
print(data.loc[~mask_utd, 'rp_joint'].mean())

