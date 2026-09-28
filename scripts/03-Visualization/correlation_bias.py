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

# Observed correlations
dfs = []
n_files = range(50,2351,50)
for n in n_files:
 print(n)
 df  = pd.read_csv(f'C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/correlations/{name}_cor_by_acc_time_{n}.csv')
 print(f'{n}: {len(df)} rows')
 dfs.append(df)

# Concatenate row-wise
data = pd.concat(dfs, axis=0, ignore_index=True)

# 24-hour rainfall
data = data[data['rain_dur_seq']==24].reset_index()

# Modelled correlations
#data = pd.read_csv(f'C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/Mullica_rain_res_vol_aorc_subset_2000.csv')
 
data_r =  pd.read_csv(f'C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/{name}_res_cor_con_rain_aorc_subset.csv')
data_ntr =  pd.read_csv(f'C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/{name}_res_cor_con_ntr_aorc_subset.csv')

#data_r_24 =  pd.read_csv(paste('C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/',name,'_res_cor_24hr_rain_ntr_con_24hr_rain_aorc_subset.csv',sep=""))
#data_ntr_24 =  pd.read_csv(paste('C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/',name,'_res_cor_24hr_rain_ntr_con_ntr_aorc_subset.csv',sep=""))
 
# Differences in correlations
data_r['bias_r_tc']       = data['cor_r_tc']       - data_r['cor_r_tc']
data_r['bias_r_non_tc']   = data['cor_r_non_tc']   - data_r['cor_r_non_tc']

data_ntr['bias_ntr_tc']  = data['cor_ntr_tc'] - data_ntr['cor_ntr_tc']
data_ntr['bias_ntr_non_tc']  = data['cor_ntr_non_tc'] - data_ntr['cor_ntr_non_tc']

plot_data = pd.DataFrame({
    'r_tc': data_r['bias_r_tc'].values,
    'r_non_tc': data_r['bias_r_non_tc'].values,
    'ntr_tc': data_ntr['bias_ntr_tc'].values,
    'ntr_non_tc': data_ntr['bias_ntr_non_tc'].values
})



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





## Map of relative errors

#Variables to plot
col_names = ['r_tc', 'ntr_tc', 'r_non_tc', 'ntr_non_tc']

# Create figure
# Figure
fig = plt.figure(figsize=(8,6), constrained_layout=False)

# Grid
gs = fig.add_gridspec(2, 3, width_ratios=[1.15, 2.0, 2.0],wspace=0.08, hspace=0.08)

#Create box plots
tc_data = [plot_data['r_tc'], plot_data['ntr_tc']]
non_tc_data = [plot_data['r_non_tc'],plot_data['ntr_non_tc']]

box_colors = ['white', 'white']

box_ax_tc = fig.add_subplot(gs[0, 0])
draw_boxplot(box_ax_tc, tc_data, box_colors, y_min=-1, y_max=1)


box_ax_non_tc = fig.add_subplot(gs[1, 0])
draw_boxplot(box_ax_non_tc,non_tc_data, box_colors, y_min=-1, y_max=1)

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
vmax = max(abs(plot_data[col_names].min().min()), abs(plot_data[col_names].max().max()))
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
        c=plot_data[col_names[i]],
        s=1,
        cmap=custom_cmap,
        vmin=vmin,
        vmax=vmax,
        transform=ccrs.PlateCarree()
    )

# Column titles
col_titles = ['Rainfall-conditioned', 'NTR-conditioned']  # or whatever your columns represent
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
plt.savefig(f'C:/Users/ro327497/OneDrive - University of Central Florida/Desktop/conferences/Coastal cluster/{name}_correlation_w_boxplots.png', dpi=300, bbox_inches='tight')

plt.show()

