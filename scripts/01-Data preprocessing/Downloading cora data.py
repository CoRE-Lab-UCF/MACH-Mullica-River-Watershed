#To run activate 'cora' in Anaconda Prompt: conda activate cora
#This was created by cloning the github repository: 
#    git clone https://github.com/NOAA-CO-OPS/CORA-Coastal-Ocean-ReAnalysis-CORA.git
#Then:
#    conda env create -f environment.yml
#Import relevent packages
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import dask
import intake
import xarray as xr
import holoviews as hv
import geoviews as gv
import hvplot.xarray
import holoviews.operation.datashader as dshade
from bokeh.models import DatetimeTickFormatter, HoverTool
import cmocean

#Setting up Python environment to render interactive visualizations
hv.extension('bokeh')

from bokeh.resources import INLINE
import bokeh.io
from bokeh import *
bokeh.io.output_notebook(INLINE)

#Acceesses a .yml file located on the National Center (NODD) which shows which CORA output files are downloadable 
catalog = intake.open_catalog("s3://noaa-nos-cora-pds/CORA_V1.1_intake.yml",storage_options={'anon':True})
list(catalog)

#to_dask() converts the dataset to a Dask DataArray (lazy loading parrel computing)
ds = catalog["CORA-V1.1-fort.63-timeseries"].to_dask()
ds

#Extract data from a specified time and create the plot

# find the indices of the points in (x,y) closest to the points in (xi,yi)
def nearxy(x,y,xi,yi):
    ind = np.ones(len(xi),dtype=int)
    for i in range(len(xi)):
        dist = np.sqrt((x-xi[i])**2+(y-yi[i])**2)
        ind[i] = dist.argmin()
    return ind

#Lat/long of point of interest
lat = 32.775
lon = -79.9239

#Get water levels of node closest to this point
ind = nearxy(ds['x'].values,ds['y'].values,[lon], [lat])

#Times
start="1989-01-01"
end="1989-12-31"

tickfmt = DatetimeTickFormatter(years="%m/%d/%Y", months="%m/%d/%Y", days = '%m/%d/%Y', hourmin = '%H:%M')
tooltips = [
    ("time", "@time{%F %T}"),
    ("water level", "@zeta"),
]
hover = HoverTool(tooltips=tooltips,formatters={
        '@time': 'datetime'})

zeta_tslice = ds['zeta'][:,ind].sel(time=slice(start, end)).compute()

plot = zeta_tslice.hvplot(x='time', grid=True, xformatter=tickfmt, tools=[hover], width=1000, height=700)
hv.save(plot, 'F:OneDrive - University of Central Florida/Documents/CONUS/Data/zeta_plot.html', backend='bokeh')

#Plot grid
import matplotlib.pyplot as plt

x = ds['x'].values
y = ds['y'].values

plt.figure(figsize=(8, 6))
plt.scatter(x, y, s=2, color='blue')  # s=2 keeps the points small and neat
plt.title('Scatter Plot of Grid Nodes')
plt.xlabel('X')
plt.ylabel('Y')
plt.grid(True)
plt.axis('equal')  # Optional: makes x and y scales equal
plt.show()

#As for AORC
#NE called NJ: latitude=slice(33, 48), longitude=slice(-90, -60)
##SE latitude=slice(20, 33), longitude=slice(-90, -60)

#Obtaining points within a grid
xmin, xmax = -90, -60
ymin, ymax = 20, 33

mask = (ds.x >= xmin) & (ds.x <= xmax) & (ds.y >= ymin) & (ds.y <= ymax)

# Apply the mask to get the node indices
# Apply the mask to get the node indices
mask_computed = mask.compute().values

# Use nonzero to get the indices where the mask is True
import numpy as np
node_indices = np.nonzero(mask_computed)[0]

# Print the indices
print(node_indices)

# Then extract zeta for all those nodes over time:
zeta_subset = ds.zeta[:, node_indices]

x = zeta_subset['x'].values
y = zeta_subset['y'].values

plt.figure(figsize=(8, 6))
plt.scatter(x, y, s=2, color='blue')  # s=2 keeps the points small and neat
plt.title('Scatter Plot of Grid Nodes')
plt.xlabel('X')
plt.ylabel('Y')
plt.grid(True)
plt.axis('equal')  # Optional: makes x and y scales equal
plt.show()


# Save with compression
file_name = f'D:/OneDrive - University of Central Florida/CORA_SE.nc'  # Using an index for unique file names
encoding = {'zeta': {"zlib": True, "complevel": 4, "shuffle": True}}
zeta_subset.to_netcdf(file_name, engine="netcdf4", format="NETCDF4", encoding=encoding)
