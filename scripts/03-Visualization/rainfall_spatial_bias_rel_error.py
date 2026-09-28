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

data['rel_bias_40_non_tc'].mean()
data['rel_bias_40_tc'].mean()


##Box plots
from matplotlib.patches import Rectangle

# ── Function: compute box stats for one vector ─────────────────
def box_stats(x):
    x = np.asarray(x)
    x = x[~np.isnan(x)]

    q1 = np.quantile(x, 0.25)
    med = np.median(x)
    q3 = np.quantile(x, 0.75)

    iqr = q3 - q1

    # Lower and upper whisker ends
    w_low = np.min(x[x >= q1 - 1.5 * iqr])
    w_high = np.max(x[x <= q3 + 1.5 * iqr])

    # Outliers
    out = x[(x < w_low) | (x > w_high)]

    return {'q1': q1, 'med': med, 'q3': q3, 'w_low': w_low, 'w_high': w_high, 'out': out}

data_list = [
    data['diff_10_tc'],
    data['diff_40_tc'],
    data['diff_10_non_tc'],
    data['diff_40_non_tc']
]

#x-axis labels
labels = ['10-year', '40-year']

#Box colors
box_cols = ['red','blue','red','blue']

# Compute stats for every dataset
stats = [box_stats(x) for x in data_list]

# ── Set up blank plot ──────────────────────────────────────────
all_vals = np.concatenate([np.asarray(x)[~np.isnan(np.asarray(x))] for x in data_list])

fig, ax = plt.subplots(figsize=(8, 6))

ax.set_xlim(0.5, 9.5)
ax.set_ylim( np.min(all_vals) - 5, np.max(all_vals) + 5)

# Horizontal zero line
ax.axhline(0, color='black', linewidth=1)

# Axis labels
ax.set_xlabel('Event', labelpad=10)
ax.set_ylabel('Bias (mm)', labelpad=10)

# Custom x-axis labels
ax.set_xticks([1, 3])
ax.set_xticklabels(labels)

# Grid
ax.grid(axis='y', linestyle=':', color='gray', alpha=0.5)


# ── Draw each box-and-whisker manually ─────────────────────────

box_width = 0.35 / 2
cap_width = 0.15 / 2

x_positions = [
    0.75, 1.25,
    2.75, 3.25]

for i in range(len(data_list)):

    s = stats[i]
    cx = x_positions[i]

    # 1. Box (Q1 to Q3)
    rect = Rectangle((cx - box_width, s['q1']), 2 * box_width, s['q3'] - s['q1'], facecolor=box_cols[i], edgecolor='black', linewidth=1.5)
    ax.add_patch(rect)

    # 2. Median line
    ax.plot([cx - box_width, cx + box_width], [s['med'], s['med']], color='black', linewidth=2.5)

    # 3. Whisker lines
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

    # 4. Whisker caps
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

    # 5. Outliers
    if len(s['out']) > 0:

        ax.scatter(
            np.full(len(s['out']), cx),
            s['out'],
            s=30,
            marker='o',
            color=box_cols[i],
            edgecolor='none'
        )


plt.tight_layout()
plt.show()


# Combined boxplots and bias 

#Variables to plot
col_names = ['diff_10_tc', 'diff_40_tc', 'diff_10_non_tc', 'diff_40_non_tc']


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
    fraction=0.025,       # was 0.035 — smaller bar
    pad=0.03,            # tighter to panels
    shrink=0.9           # shortens the bar vertically
    )

cbar.set_label('Bias (mm)', fontsize=11)

# save
plt.savefig(f'C:/Users/ro327497/OneDrive - University of Central Florida/Desktop/conferences/Coastal cluster/{name}_rainfall_spatial_bias.png', dpi=300, bbox_inches=None)

plt.show()








## Relative errors

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
    fraction=0.025,       # was 0.035 — smaller bar
    pad=0.03,            # tighter to panels
    shrink=0.9           # shortens the bar vertically
    )

cbar.set_label('Relative bias (%)', fontsize=11)

# save
plt.savefig(f'C:/Users/ro327497/OneDrive - University of Central Florida/Desktop/conferences/Coastal cluster/{name}_rainfall_spatial_rel_error_v1.png', dpi=300, bbox_inches=None)

plt.show()