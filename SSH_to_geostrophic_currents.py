import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import math
from mpl_toolkits.basemap import Basemap

rr = Dataset('ssh_20071127.nc', 'r', format='NETCDF3')
ssh = rr.variables['ssh'][:]
lat = rr.variables['latitude'][:]
lon = rr.variables['longitude'][:]
L=40075000

u=np.zeros((len(lat)-1,len(lon)-1)); v=np.zeros((len(lat)-1,len(lon)-1))
f=[2*7.2921*1e-05*np.sin(math.radians(lat[i])) for i in range(len(lat))]
p=np.array([9.81*(ssh[i,j]-0)+1013.25 for i in range(len(ssh)) for j in range(len(ssh.T))]).reshape(len(lat),len(lon))
v=np.array([(p[i,j+1]-p[i,j])/((lon[j+1]-lon[j])*((1/360.)*L*np.cos(math.radians((lat[i+1]+lat[i])*0.5)))*f[i]) for i in range(len(lat)-1) for j in range(len(lon)-1)]).reshape(len(lat)-1,len(lon)-1)
u=np.array([(p[i+1,j]-p[i,j])/((lat[i+1]-lat[i])*111200*f[i]) for i in range(len(lat)-1) for j in range(len(lon)-1)]).reshape(len(lat)-1,len(lon)-1)
w=(u**2+v**2)**0.5
plt.imshow(np.flipud(w))
plt.show()
