from pylab import *; from netCDF4 import Dataset

var='v' # 'temp'
lat=80
lon=75
nc = Dataset('ocean_nest_calvi_his_surf.nc', 'r', format='netCDF4')
temp=nc.variables[var][:]

# COmputing temporal averages for a variable

av=np.zeros((len(temp[0,0,:,0]),len(temp[0,0,0,:])))
for i in range(len(av)):
    for j in range(len(av.T)):
        av[i,j]=mean(temp[:,0,i,j])

# Computing three-dimensional martrix for a variable (unique value minus a corresponding average value)
cov=np.zeros((len(temp),len(av),len(av.T)))
for i in range(len(temp)):
    for j in range(len(temp[0,0,:,0])):
        for k in range(len(temp[0,0,0,:])):
            cov[i,j,k]=temp[i,0,j,k]-av[j,k]


# Computing the 3D covariance matrix for the chosen grid point
covv=np.zeros((len(temp),len(av),len(av.T)))
cov_u=np.zeros((len(temp)))
for i in range(len(temp)):
	#covv[i]=(cov[i].flatten()*cov[i,lat,lon]).reshape(shape(temp)[-2],shape(temp)[-1])
	cov_u[i]=nc.variables['u'][i,0,lat,lon]-mean(nc.variables['u'][:,0,lat,lon])
	covv[i]=(cov[i].flatten()*cov_u[i]).reshape(shape(temp)[-2],shape(temp)[-1])

"""
# Computing the 2D covariance matrix for each spatial grid point
covar=np.zeros((len(av),len(av.T)))
for i in range(len(covar)):
    for j in range(len(covar.T)):
        covar[i,j]=mean(covv[:,i,j])
"""
# Computing the standard deviation for each spatial grid point
stdev=np.zeros((len(av),len(av.T)))
#stdfix=std(temp[:,0,lat,lon])
stdfix=std(nc.variables['u'][:,0,lat,lon])
for i in range(len(av)):
    for j in range(len(av.T)):
        stdev[i,j]=mean(covv[:,i,j])/(std(temp[:,0,i,j])*stdfix)

stdev_mask=ma.masked_where(np.isnan(stdev),stdev)
a=plt.pcolormesh(stdev_mask)
colorbar()
plot(lon,lat,'o')

show()
