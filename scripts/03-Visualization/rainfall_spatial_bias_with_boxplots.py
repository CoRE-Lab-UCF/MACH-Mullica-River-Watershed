# -*- coding: utf-8 -*-
"""
Created on Tue Aug 18 15:45:34 2026

@author: ro327497
"""

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


##Loop through return levels
#Correltations
data = pd.read_csv(f'C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/Mullica_rain_res_vol_aorc_subset_2000.csv')
         
#
data['diff_10_tc'] = data['obs_rp10_tc'] - data['rp10_50_tc']
data['diff_40_tc'] = data['obs_rp40_tc'] - data['rp40_50_tc']
data['rel_bias_10_tc'] = 100*(data['obs_gpd_rp10_tc'] - data['gpd_rp10_50_tc'])/ data['obs_gpd_rp10_tc']
data['rel_bias_40_tc'] = 100*(data['obs_gpd_rp40_tc'] - data['gpd_rp40_50_tc'])/ data['obs_gpd_rp40_tc']
#data['diff_50_tc'] = data['obs_gpd_rp50_tc'] - data['gpd_rp50_50_tc']
#data['diff_100_tc'] = data['obs_gpd_rp100_tc'] - data['gpd_rp100_50_tc'] 

data['diff_10_non_tc'] = data['obs_rp10_non_tc'] - data['rp10_50_non_tc']
data['diff_40_non_tc'] = data['obs_rp40_non_tc'] - data['rp40_50_non_tc']
data['rel_bias_10_non_tc'] = 100*(data['obs_gpd_rp10_non_tc'] - data['gpd_rp10_50_non_tc']) /  data['obs_gpd_rp10_non_tc']
data['rel_bias_40_non_tc'] = 100*(data['obs_gpd_rp40_non_tc'] - data['gpd_rp40_50_non_tc']) /  data['obs_gpd_rp40_non_tc']
#data['diff_100_non_tc'] = data['obs_gpd_rp100_non_tc'] - data['gpd_rp100_50_non_tc'] 

#Min and max of data
data[['diff_10_tc', 'diff_40_tc', 'rel_bias_10_tc', 'rel_bias_40_tc']].max() 
data[['diff_10_non_tc', 'diff_40_non_tc', 'rel_bias_10_non_tc', 'rel_bias_40_non_tc']].max() 

data[['diff_10_tc', 'diff_40_tc', 'rel_bias_10_tc', 'rel_bias_40_tc']].min() 
data[['diff_10_non_tc', 'diff_40_non_tc', 'rel_bias_10_non_tc', 'rel_bias_40_non_tc']].min() 

data[['rel_bias_40_tc', 'rel_bias_40_non_tc']].mean()
  

 
##Box plots
from matplotlib.patches import Rectangle

def box_stats(x):
    x = np.asarray(x)
    x = x[~np.isnan(x)]

    q1 = np.quantile(x, 0.25)
    med = np.median(x)
    q3 = np.quantile(x, 0.75)

    iqr = q3 - q1

    w_low = np.min(x[x >= q1 - 1.5 * iqr])
    w_high = np.max(x[x <= q3 + 1.5 * iqr])

    out = x[(x < w_low) | (x > w_high)]

    return {
        'q1': q1,
        'med': med,
        'q3': q3,
        'w_low': w_low,
        'w_high': w_high,
        'out': out
    }


def draw_boxplot(ax, data_list, colors, y_min, y_max):

    stats = [box_stats(x) for x in data_list]

    box_width = 0.25
    cap_width = 0.12

    x_positions = [1, 2]

    all_vals = np.concatenate([
        np.asarray(x)[~np.isnan(np.asarray(x))]
        for x in data_list
    ])

    for i, s in enumerate(stats):

        cx = x_positions[i]

        # Box
        rect = Rectangle(
            (cx - box_width, s['q1']),
            2 * box_width,
            s['q3'] - s['q1'],
            facecolor=colors[i],
            edgecolor='black',
            linewidth=1.5
        )

        ax.add_patch(rect)

        # Median
        ax.plot(
            [cx - box_width, cx + box_width],
            [s['med'], s['med']],
            color='black',
            linewidth=2
        )

        # Whiskers
        ax.plot(
            [cx, cx],
            [s['q1'], s['w_low']],
            color='black',
            linewidth=1.5
        )

        ax.plot(
            [cx, cx],
            [s['q3'], s['w_high']],
            color='black',
            linewidth=1.5
        )

        # Caps
        ax.plot(
            [cx - cap_width, cx + cap_width],
            [s['w_low'], s['w_low']],
            color='black',
            linewidth=1.5
        )

        ax.plot(
            [cx - cap_width, cx + cap_width],
            [s['w_high'], s['w_high']],
            color='black',
            linewidth=1.5
        )

        # Outliers
        if len(s['out']) > 0:

            ax.scatter(
                np.full(len(s['out']), cx),
                s['out'],
                s=20,
                color=colors[i],
                edgecolor='none'
            )
    
    ax.set_ylim(y_min, y_max)
    ax.set_xlim(0.5, 2.5)

    #ax.set_xticks([1, 2])
    #ax.set_xticklabels(['10-year', '40-year'])

    ax.axhline(
        0,
        color='black',
        linewidth=1
    )

    ax.grid(
        axis='y',
        linestyle=':',
        color='gray',
        alpha=0.5
    )

    ax.set_ylabel('Bias (mm)', fontsize=10, labelpad=1)



# Combined boxplots and bias 

#Variables to plot
col_names = ['diff_10_tc', 'diff_40_tc', 'diff_10_non_tc', 'diff_40_non_tc']


# Create figure
# Figure
fig = plt.figure(figsize=(8,6), constrained_layout=False)

# Grid
gs = fig.add_gridspec(2, 3, width_ratios=[1.15, 2.0, 2.0],wspace=0.08, hspace=0.08)

#Create box plots
tc_data = [data['diff_10_tc'], data['diff_40_tc']]
non_tc_data = [data['diff_10_non_tc'],data['diff_40_non_tc']]

box_colors = ['white', 'white']

box_ax_tc = fig.add_subplot(gs[0, 0])
draw_boxplot(box_ax_tc, tc_data, box_colors, y_min=-30, y_max=82)


box_ax_non_tc = fig.add_subplot(gs[1, 0])
draw_boxplot(box_ax_non_tc,non_tc_data, box_colors, y_min=-30, y_max=82)

for ax in [box_ax_tc, box_ax_non_tc]:
    pos = ax.get_position()

    new_height = pos.height * 0.75
    new_y = pos.y0 + (pos.height - new_height) / 2

    ax.set_position([
        pos.x0,
        new_y,
        pos.width,
        new_height
    ])
 
box_ax_tc.set_xticklabels([])   
box_ax_tc.tick_params(axis='x', length=0)  
box_ax_non_tc.tick_params(axis='x', length=0)  
box_ax_non_tc.set_xticks([1, 2])
box_ax_non_tc.set_xticklabels(['10-year', '40-year'])
    
# store axes as 2D
axes = [[None for _ in range(2)] for _ in range(2)]

# custom colormap
# Create a symetric custom colormap
custom_cmap = plt.get_cmap('RdBu_r')
vmax = max(abs(data[col_names].min().min()), abs(data[col_names].max().max()))
vmin = -vmax

# Create a ScalarMappable for the colorbar
norm = plt.Normalize(vmin=vmin, vmax=vmax)
sm = plt.cm.ScalarMappable(cmap=custom_cmap, norm=norm)
sm.set_array([])

#Define terrain
stamen_terrain = cimgt.StadiaMapsTiles(apikey='1dabf9e0-8140-4bfc-9327-2c1ef91e4557', style='stamen_terrain_background', resolution='@2x')

#Compute min and max of color bar
# -------------------------
# PLOTTING LOOP
# -------------------------
for i in range(4):

    row = i // 2
    col = i % 2

    ax = fig.add_subplot(gs[row, col+1], projection=stamen_terrain.crs)
    axes[row][col] = ax

    # extent
    ax.set_extent([-75, -74, 39, 40], crs=ccrs.PlateCarree())

    # add basemap
    ax.add_image(stamen_terrain, 8)
    ax.add_feature(cfeature.COASTLINE, linewidth=0.5)
    ax.add_feature(cfeature.STATES, linewidth=0.5)
    ax.add_feature(cfeature.RIVERS)
    ax.add_feature(cfeature.LAKES)

    # scatter
    ax.scatter(
        aorc_subset.lon,
        aorc_subset.lat,
        c=data[col_names[i]],
        s=1,
        cmap=custom_cmap,
        vmin=vmin,
        vmax=vmax,
        transform=ccrs.PlateCarree()
    )

# Column titles
col_titles = ['10-year', '40-year']  # or whatever your columns represent
for c in range(2):
    axes[0][c].set_title(col_titles[c], fontsize=11,pad=5)

# Row titles
row_titles = ['TC', 'Non-TC']

for r in range(2):
    # get the left axis of each row
    ax = axes[r][0]
    ax.text(
        -1.15, 0.5,          # just outside left edge
        row_titles[r],
        transform=ax.transAxes,
        fontsize=11,
        va='center',
        ha='center',
        rotation=0
    )

row_sm = []

# Figure placement

plt.subplots_adjust(
    left=0.13,
    right=0.88,
    top=0.98,
    bottom=0.02
)

all_axes = [ax for row in axes for ax in row]

cbar = fig.colorbar(
    sm,
    ax=all_axes,
    orientation='vertical',
    fraction=0.035,       # was 0.035 — smaller bar
    pad=0.03,            # tighter to panels
    shrink=0.9           # shortens the bar vertically
    )

cbar.set_label('Bias (mm)', fontsize=11)

# save
plt.savefig(f'C:/Users/ro327497/OneDrive - University of Central Florida/Desktop/conferences/Coastal cluster/{name}_rainfall_spatial_bias_w_boxplots.png', dpi=300, bbox_inches='tight')

plt.show()




## Map of relative errors

#Variables to plot
col_names = ['rel_bias_10_tc', 'rel_bias_40_tc', 'rel_bias_10_non_tc', 'rel_bias_40_non_tc']

# Create figure
# Figure
fig = plt.figure(figsize=(6,6), constrained_layout=False)

# Grid
gs = fig.add_gridspec(2, 2, wspace=0.03, hspace=0.01)

# store axes as 2D
axes = [[None for _ in range(2)] for _ in range(2)]

# custom colormap
# Create a symetric custom colormap
custom_cmap = plt.get_cmap('RdBu_r')
vmax = max(abs(data[col_names].min().min()), abs(data[col_names].max().max()))
vmin = -vmax

# Create a ScalarMappable for the colorbar
norm = plt.Normalize(vmin=vmin, vmax=vmax)
sm = plt.cm.ScalarMappable(cmap=custom_cmap, norm=norm)
sm.set_array([])

#Define terrain
stamen_terrain = cimgt.StadiaMapsTiles(apikey='1dabf9e0-8140-4bfc-9327-2c1ef91e4557', style='stamen_terrain_background', resolution='@2x')

#Compute min and max of color bar
# -------------------------
# PLOTTING LOOP
# -------------------------
for i in range(4):

    row = i // 2
    col = i % 2

    ax = fig.add_subplot(gs[row, col], projection=stamen_terrain.crs)
    axes[row][col] = ax

    # extent
    ax.set_extent([-75, -74, 39, 40], crs=ccrs.PlateCarree())

    # add basemap
    ax.add_image(stamen_terrain, 8)
    ax.add_feature(cfeature.COASTLINE, linewidth=0.5)
    ax.add_feature(cfeature.STATES, linewidth=0.5)
    ax.add_feature(cfeature.RIVERS)
    ax.add_feature(cfeature.LAKES)

    # scatter
    ax.scatter(
        aorc_subset.lon,
        aorc_subset.lat,
        c=data[col_names[i]],
        s=1,
        cmap=custom_cmap,
        vmin=vmin,
        vmax=vmax,
        transform=ccrs.PlateCarree()
    )

# Column titles
col_titles = ['10-year', '40-year']  # or whatever your columns represent
for c in range(2):
    axes[0][c].set_title(col_titles[c], fontsize=11,pad=3)

# Row titles
row_titles = ['TC', 'Non-TC']

for r in range(2):
    # get the left axis of each row
    ax = axes[r][0]
    ax.text(
        -0.2, 0.5,          # just outside left edge
        row_titles[r],
        transform=ax.transAxes,
        fontsize=11,
        va='center',
        ha='center',
        rotation=0
    )

row_sm = []

# Figure placement

plt.subplots_adjust(
    left=0.13,
    right=0.88,
    top=0.98,
    bottom=0.02
)

all_axes = [ax for row in axes for ax in row]

cbar = fig.colorbar(
    sm,
    ax=all_axes,
    orientation='vertical',
    fraction=0.035,       # was 0.035 — smaller bar
    pad=0.03,            # tighter to panels
    shrink=0.9           # shortens the bar vertically
    )

cbar.set_label('Relative bias (%)', fontsize=11)

# save
plt.savefig(f'C:/Users/ro327497/OneDrive - University of Central Florida/Desktop/conferences/Coastal cluster/{name}_rainfall_spatial_rel_error.png', dpi=300, bbox_inches=None)

plt.show()




## Map of 10 and 40 year rainfall totals

#Variables to plot
col_names = ['obs_rp10_tc', 'obs_rp40_tc', 'obs_rp10_non_tc', 'obs_rp40_non_tc']

# Create figure
# Figure
fig = plt.figure(figsize=(6,6), constrained_layout=False)

# Grid
gs = fig.add_gridspec(2, 2, wspace=0.03, hspace=0.01)

# store axes as 2D
axes = [[None for _ in range(2)] for _ in range(2)]

# custom colormap
# Create a symetric custom colormap
custom_cmap = plt.get_cmap('GnBu')
vmax = data[col_names].max().max()
vmin = 0

# Create a ScalarMappable for the colorbar
norm = plt.Normalize(vmin=vmin, vmax=vmax)
sm = plt.cm.ScalarMappable(cmap=custom_cmap, norm=norm)
sm.set_array([])

#Define terrain
stamen_terrain = cimgt.StadiaMapsTiles(apikey='1dabf9e0-8140-4bfc-9327-2c1ef91e4557', style='stamen_terrain_background', resolution='@2x')

#Compute min and max of color bar
# -------------------------
# PLOTTING LOOP
# -------------------------
for i in range(4):

    row = i // 2
    col = i % 2

    ax = fig.add_subplot(gs[row, col], projection=stamen_terrain.crs)
    axes[row][col] = ax

    # extent
    ax.set_extent([-75, -74, 39, 40], crs=ccrs.PlateCarree())

    # add basemap
    ax.add_image(stamen_terrain, 8)
    ax.add_feature(cfeature.COASTLINE, linewidth=0.5)
    ax.add_feature(cfeature.STATES, linewidth=0.5)
    ax.add_feature(cfeature.RIVERS)
    ax.add_feature(cfeature.LAKES)

    # scatter
    ax.scatter(
        aorc_subset.lon,
        aorc_subset.lat,
        c=data[col_names[i]],
        s=1,
        cmap=custom_cmap,
        vmin=vmin,
        vmax=vmax,
        transform=ccrs.PlateCarree()
    )

# Column titles
col_titles = ['10-year', '40-year']  # or whatever your columns represent
for c in range(2):
    axes[0][c].set_title(col_titles[c], fontsize=11,pad=3)

# Row titles
row_titles = ['TC', 'Non-TC']

for r in range(2):
    # get the left axis of each row
    ax = axes[r][0]
    ax.text(
        -0.2, 0.5,          # just outside left edge
        row_titles[r],
        transform=ax.transAxes,
        fontsize=11,
        va='center',
        ha='center',
        rotation=0
    )

row_sm = []

# Figure placement

plt.subplots_adjust(
    left=0.13,
    right=0.88,
    top=0.98,
    bottom=0.02
)

all_axes = [ax for row in axes for ax in row]

cbar = fig.colorbar(
    sm,
    ax=all_axes,
    orientation='vertical',
    fraction=0.035,       # was 0.035 — smaller bar
    pad=0.03,            # tighter to panels
    shrink=0.9           # shortens the bar vertically
    )


cbar.set_label('Total rainfall (mm)', fontsize=11)

# save
plt.savefig(f'C:/Users/ro327497/OneDrive - University of Central Florida/Desktop/conferences/Coastal cluster/{name}_rainfall_spatial_rainfall_totals.png', dpi=300, bbox_inches=None)

plt.show()