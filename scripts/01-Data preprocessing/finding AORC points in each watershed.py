#Install geopandas
#conda install geopandas (done)

#Import relevent packages
import geopandas as gpd

#Load the shape file
watershed_10 = gpd.read_file("F:/OneDrive - University of Central Florida/Documents/USACE - Will/Maps/WBD_03_HU2_Shape/Shape/WBDHU10.shp")

#Exploring the .shp file
print(watershed_10.columns)      # Check all available fields
print(watershed_10.head())       # Preview the first few rows
watershed_10.plot()              # Optional: quick visual plot

print(watershed_10.states) 

#Watersheds in NJ are 
nj_watersheds = watershed_10[watershed_10["states"] == "NJ"]

#Number of watersheds
nj_watersheds.shape

#Finding which AORC points are in NJ

# Find which polygons contain this point
containing_polygons = []

for i, polygon in enumerate(states_geometries):
    if polygon.contains(point):
        containing_polygons.append((i, polygon))

# Print the results
if containing_polygons:
    print(f"Found {len(containing_polygons)} polygons containing point ({point_x}, {point_y}):")
    for i, polygon in containing_polygons:
        print(f"  Polygon index: {i}")
        print(f"  Polygon type: {type(polygon)}")
        print(f"  Polygon area: {polygon.area:.6f} square degrees")
        #print(f"  Polygon coordinates count: {len(polygon.exterior.coords)}")
        #print(f"  First coordinate: {list(polygon.exterior.coords)[0]}")
        print("-" * 40)
    
    # Visualize the containing polygons
    fig = plt.figure(figsize=(10, 6))
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
    
    # Add all state boundaries in light gray
    ax.add_feature(states, edgecolor='gray', facecolor='none', linewidth=0.5)
    
    # Highlight the containing polygons
    for i, polygon in containing_polygons:
        ax.add_geometries([polygon], crs=ccrs.PlateCarree(), 
                         facecolor='red', edgecolor='black', alpha=0.7)
    
    # Mark the point
    ax.plot(point_x, point_y, 'bo', markersize=8, transform=ccrs.PlateCarree())
    
    # Set the map extent to focus on the area around the point
    padding = 10  # degrees
    ax.set_extent([point_x-padding, point_x+padding, 
                  point_y-padding, point_y+padding], ccrs.PlateCarree())
    
    plt.title(f"Polygons containing point ({point_x}, {point_y})")
    plt.show()
else:
    print(f"No polygons found containing point ({point_x}, {point_y})")
    
    
#Extract polygon
coordinates = [list(polygon.exterior.coords) for polygon in states_geometries[80].geoms]

print("Coordinates:", coordinates)
# Reproject watershed polygons to match points CRS
nj_watersheds = nj_watersheds.to_crs(points_gdf.crs)

# Perform spatial join
points_in_watershed = gpd.sjoin(points_gdf, nj_watersheds, how='inner', predicate='within')

# Get the indices of points that are inside the watershed
points_indices_in_watershed = points_in_watershed.index


####

#Read in data
import xarray as xr
import geopandas as gpd
from shapely.geometry import Point
import pandas as pd
import numpy as np

# === Step 1: Load the NetCDF file ===
ds = xr.open_dataset("F:/OneDrive - University of Central Florida/Documents/CONUS/Data/AORC/AORC_1km_NE_1979.nc")

# Replace with your actual variable names for latitude and longitude
lat = ds['latitude'].values
lon = ds['longitude'].values

# If it's a meshgrid format, flatten for point-wise access
lon2d, lat2d = np.meshgrid(lon, lat)
lat_flat = lat2d.flatten()
lon_flat = lon2d.flatten()

# === Step 2: Convert to GeoDataFrame of points ===
points_df = pd.DataFrame({'lat': lat_flat, 'lon': lon_flat})
geometry = [Point(xy) for xy in zip(points_df.lon, points_df.lat)]
points_gdf = gpd.GeoDataFrame(points_df, geometry=geometry, crs="EPSG:4326")

# === Step 3: Load watershed shapefile ===
watersheds = gpd.read_file('F:/OneDrive - University of Central Florida/Documents/WBD/WBD_National_GDB.gdb', layer='WBDHU8')

# Ensure same CRS
watersheds = watersheds.to_crs(points_gdf.crs)

# === Step 4: Spatial join — assign each point to a watershed ===
# This will only keep points that fall inside any watershed polygon
joined = gpd.sjoin(points_gdf, watersheds, how="inner", predicate="within")

# joined now contains columns from both points and matching watershed polygons
# For example, if your watershed file had a "Name" field:
assigned_points = joined[['lat', 'lon', 'tnmid', 'states' , 'name']]

#Watersheds in NJ are 
nj_points = assigned_points[assigned_points['states'].str.contains('NJ')]

#Names of watersheds at least partially in NJ
nj_points['name'].unique()

#Number of watershedsat least partially in NJ (15)
len(nj_points['name'].unique())

#AORC.nc dataframe index of points in nj_points
matched = gpd.sjoin(joined, points_gdf.reset_index(), how='left', predicate='within')

#Add as column to nj_points GeoDataFrame
nj_points['aorc_id'] = matched['index']


import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import cartopy.io.img_tiles as cimgt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import xarray as xr
import numpy as np

#Download map
stamen_terrain = cimgt.StadiaMapsTiles(apikey='1dabf9e0-8140-4bfc-9327-2c1ef91e4557', style='stamen_terrain_background', resolution='@2x')

#Single map
fig = plt.figure(figsize=(10, 8))

cmap = plt.get_cmap('plasma_r')

ax = fig.add_subplot(projection=stamen_terrain.crs)
# set miami-dade county as map extent
ax.set_extent([-73, -75.9157, 38.9, 41.4], crs=ccrs.Geodetic())
# Add the Stamen data at zoom level 8.
ax.add_image(stamen_terrain, 8)
ax.add_feature(cfeature.COASTLINE)
ax.add_feature(cfeature.STATES)
ax.add_feature(cfeature.RIVERS)
ax.add_feature(cfeature.LAKES)
ax.add_feature(cfeature.OCEAN)
ax.add_feature(cfeature.LAND)

import matplotlib.colors as mcolors

watershed_names = nj_points['name'].unique()
cmap = plt.cm.get_cmap('tab20', len(watershed_names))
colors = {name: cmap(i) for i, name in enumerate(watershed_names)}

for name in watershed_names:
    subset = nj_points[nj_points['name'] == name]
    ax.scatter(
       subset.lon, 
       subset.lat, 
       color=colors[name], 
       s=10, 
       label=name,
       transform=ccrs.PlateCarree()  # important for cartopy plots
   )
ax.legend(loc='lower right', fontsize='small')

plt.show()


####Cora in each watershed

# === Step 1: Load the NetCDF file ===
ds = xr.open_dataset("F:/OneDrive - University of Central Florida/Documents/CONUS/Data/CORA/CORA_NJ.nc")

# Replace with your actual variable names for latitude and longitude (not a mesh)
lon = ds['x'].values
lat = ds['y'].values

# === Step 2: Convert to GeoDataFrame of points ===
points_df = pd.DataFrame({'lat': lat, 'lon': lon})
geometry = [Point(xy) for xy in zip(points_df.lon, points_df.lat)]
points_gdf = gpd.GeoDataFrame(points_df, geometry=geometry, crs="EPSG:4326")

# === Step 3: Load watershed shapefile ===
watersheds = gpd.read_file('F:/OneDrive - University of Central Florida/Documents/WBD/WBD_National_GDB.gdb', layer='WBDHU8')

# Ensure same CRS
watersheds = watersheds.to_crs(points_gdf.crs)

# === Step 4: Spatial join — assign each point to a watershed ===
# This will only keep points that fall inside any watershed polygon
joined = gpd.sjoin(points_gdf, watersheds, how="inner", predicate="within")

# joined now contains columns from both points and matching watershed polygons
# For example, if your watershed file had a "Name" field:
assigned_points = joined[['lat', 'lon', 'tnmid', 'states' , 'name']]

print(assigned_points.head())

#Watersheds in NJ are 
watersheds_nj = watersheds[watersheds['states'].str.contains('NJ')]

#Centroids of NJ watersheds
# Reproject to a suitable projected CRS (example: EPSG:32618 for UTM zone 18N)
watersheds_nj_proj = watersheds_nj.to_crs(epsg=32618)

# Calculate centroids
watersheds_nj_proj['centroid'] = watersheds_nj_proj.geometry.centroid

# If needed, convert centroids back to lat/lon
watersheds_nj['centroid_latlon'] = watersheds_nj_proj['centroid'].to_crs(epsg=4326)

#Find cora node with closest location to centroids
points_proj = points_gdf.to_crs(epsg=32618) #Project ponts to get distance in meters
watersheds_nj_proj['centroid']

# For each target point, find the nearest candidate using Euclidean distance
nearest_points_idx = []
nearest_points = []
distances = []

for target_geom in watersheds_nj_proj['centroid']:
    dists = points_proj.geometry.distance(target_geom)
    min_idx = dists.idxmin()
    nearest_points_idx.append(min_idx)
    nearest_points.append(points_proj.geometry[min_idx])
    distances.append(dists[min_idx])

# Add to target_points GeoDataFrame
watersheds_nj_proj['nearest_point_idx'] = nearest_points_idx
watersheds_nj_proj['nearest_point'] = nearest_points
watersheds_nj_proj['distance'] = distances

# Optional: convert back to WGS84 for mapping/export
watersheds_final = watersheds_nj_proj.to_crs(epsg=4326)

print(watersheds_final['nearest_point_idx'])


#Download map
stamen_terrain = cimgt.StadiaMapsTiles(apikey='1dabf9e0-8140-4bfc-9327-2c1ef91e4557', style='stamen_terrain_background', resolution='@2x')

#Single map
fig = plt.figure(figsize=(10, 8))

cmap = plt.get_cmap('plasma_r')

ax = fig.add_subplot(projection=stamen_terrain.crs)
# set miami-dade county as map extent
ax.set_extent([-73, -75.9157, 38.9, 41.4], crs=ccrs.Geodetic())
# Add the Stamen data at zoom level 8.
ax.add_image(stamen_terrain, 8)
ax.add_feature(cfeature.COASTLINE)
ax.add_feature(cfeature.STATES)
ax.add_feature(cfeature.RIVERS)
ax.add_feature(cfeature.LAKES)
ax.add_feature(cfeature.OCEAN)
ax.add_feature(cfeature.LAND)

import matplotlib.colors as mcolors

watershed_names = nj_points['name'].unique()
cmap = plt.cm.get_cmap('tab20', len(watershed_names))
colors = {name: cmap(i) for i, name in enumerate(watershed_names)}

for name in watershed_names:
    subset = nj_points[nj_points['name'] == name]
    ax.scatter(
       subset.lon, 
       subset.lat, 
       color=colors[name], 
       s=10, 
       label=name,
       transform=ccrs.PlateCarree()  # important for cartopy plots
   )
ax.legend(loc='lower right', fontsize='small')

#Plot CORA point cloest to centroid of each watershed
# Extract nearest points based on indices from watersheds_final
nearest_points_geom = points_gdf.geometry.iloc[watersheds_final['nearest_point_idx']]

# Use GeoPandas to plot these points (with red color and size 20)
nearest_points_geom.plot(ax=ax, color='red', markersize=20, label='Nearest CORA Points',  transform=ccrs.Geodetic())


plt.show()


#Add in CORA node index for each AORC_NJ point
#Create pyhton dictonary
ws_to_cora_idx = watersheds_final.set_index('name')['nearest_point_idx'].to_dict()

# Map AORC index to nj_points
nj_points['cora_id'] = nj_points['name'].map(ws_to_cora_idx)
