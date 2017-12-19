from pylab import *; from netCDF4 import Dataset; from numpy.linalg import inv

# The goal is to compute different drigter trajectories and see, how the flow field will very

nc = Dataset('ocean_his_surf.nc', 'r', format='netCDF4')
uu=nc.variables['u'][0,0,:,:]
vv=nc.variables['v'][0,0,:,:]
lat=nc.variables['lat_rho'][1:-1,1:-1]
lon=nc.variables['lon_rho'][1:-1,1:-1]
pn,pm=nc.variables['pn'][:],nc.variables['pm'][:]

# put them on the cetner of the grid cells

u = (uu[1:-1,0:-1]+uu[1:-1,1:])/2.
v = (vv[1:,1:-1]+vv[0:-1,1:-1])/2.

# compute the divergence
#div = (uu[1:-1,1:]-uu[1:-1,0:-1])*pm[1:-1,1:-1] + (vv[1:,1:-1]-vv[0:-1,1:-1])*pn[1:-1,1:-1]; pcolor(lon,lat,div, cmap="RdBu_r", vmin=-1*div.max()/3.,vmax=div.max()/3.);colorbar()
#ylim(lat.min(),lat.max()); xlim(lon.min(),lon.max()); show()

def pt():
	quiver(lon[::3,::3],lat[::3,::3],u[::3,::3],v[::3,::3])
	plot(x_pos[0,:],x_pos[1,:], linewidth = 3, c = 'r', label='true')
	#plot(x_noise[0,:],x_pos[1,:], label='true + noise')
	#plot(x_free[0,:],x_free[1,:], label='no correction')
	for i in range(Nens):
		plot(x_ens[0,i,:],x_ens[1,i,:])
	ylim(lat.min(),lat.max()); xlim(lon.min(),lon.max())
	#scatter(x_pos[0,:],x_pos[1,:],s=155,c = P[0,0,:]+P[1,1,:], cmap='rainbow',edgecolors = None)
	#colorbar()
	#legend(loc=1)
	show()


def crd(lon_t,lat_t):
	fract_i = (lon_t-lon[0,0])/(lon[0,1]-lon[0,0])
	fract_j = (lat_t-lat[0,0])/(lat[1,0]-lat[0,0])
	i,j = floor(fract_i),floor(fract_j)
	di,dj = fract_i-i,fract_j-j
	return i,j,di,dj

def vel(t, x):
	if x[0]>lon.max() or x[1]>lat.max() or x[0]<lon.min() or x[1]<lat.min():
		return array([0,0])
	j,i,dj,di = crd(x[0],x[1])
	u00,u10,u01,u11 = u[i,j],u[i+1,j],u[i,j+1],u[i+1,j+1]
	ub = u00 + di * (u10-u00) 
	ut = u10 + di * (u11-u10) 
	ui = ub + dj * (ut - ub)
	v00,v10,v01,v11 = v[i,j],v[i+1,j],v[i,j+1],v[i+1,j+1]
	vb = v00 + di * (v10-v00) 
	vt = v10 + di * (v11-v10) 
	vi = vb + dj * (vt - vb)
	ui = ui / np.cos(np.deg2rad(x[1])) / 1.11e5
	vi = vi / 1.11e5
	#print ui,vi
	return array([ui,vi])
	
# define runge-kutta function	
def rgk(x, t, dt):
	xp = x + dt * vel(t,x) / 2. 
	x = x + dt * vel(t,xp)
	return x

# derive tangent-linear model for rgk (by finit differencies)
def tanmod_vector(x, dx, t, dt):
	vector = ( rgk(x+dx*eps, t, dt) - rgk(x-dx*eps, t, dt) ) * 0.5 / eps	# how the perturbation dx propagates the model
	return vector

# Than more the flow DIVERGE, than larger the error
def tanmod_matrix(x, P, t, dt):
	MP = np.zeros((2,2))
	for i in range(2):
		MP[i,:] = tanmod_vector(x, P[i,:], t, dt)		# product between the tangent linear model and the covariacne matrix P
	return MP


	

# initialization

x = [8.2 , 42.5]	# lon, lat of initial point
dt = 3600 	# s
n_it = 600 	# number of iterations
eps = 1e-4
err=np.zeros((2))
R = np.eye(2)/1000. # uncertainty in the observations


err[0] = np.sqrt(R[0,0]) * np.random.randn() # error of longitude
err[1] = np.sqrt(R[1,1]) * np.random.randn() # error of latitude

# declaration of array for a priori quantities
x_pos = np.zeros((2,n_it)) # true
x_noise = np.zeros((2,n_it)) # true+noise
x_free = np.zeros((2,n_it)) # free



P = np.zeros((2,2,n_it))
K = np.zeros((2,2,n_it))

x_pos[:,0] = x
x_free[:,0] = [8,42.5]
x_noise[:,0] = [8,42.5]
P_init = np.eye(2)/100.
P[:,:,0] = P_init				# error covariance

n_obs = 20 # assimilation of each 20 obs

H = np.eye(2)	# observation operator, extracts that we observe in the system	

################
# For ensemble #
Nens = 10

Pert = np.zeros((2,Nens))
for i in range(Nens):
	Pert[:,i] = (np.matrix(x).T + np.matrix(np.sqrt(P_init))*np.matrix(np.random.randn(2)).T).flatten()

x_ens = np.zeros((2, Nens,n_it))
x_ens[:,:,0] = Pert
################

# if column of S, wheer S is a square root matrix of P


def ETKF(x_ens,H,yo,R):
	# Input: x_ens[:,:,n] = ensemble, H = operator, yo = obs, R = error covariance
	S = np.zeros((2,Nens))
	xmean = np.mean(x_ens[:,:],1) # ensemble mean
	for i in range(Nens):
		S[:,i] = (x_ens[:,i] - xmean) / np.sqrt(Nens-1)
	SS = np.matrix(S)*np.matrix(S).T # covariance of our initial position, than more members, than closer we get to our initial position
	HS = np.matrix(H)*np.matrix(S) # scale of ensemble perturbation at thelocation of the observation
	a = HS.T * inv(R) * HS # the ratio of the model error versus the observation error, generalization of the concep of signal to noise ratio
	# now we compute eigenvalues and eigenvectors
	lamda = np.linalg.eigh(a)[0] # eigenvalues
	U = np.linalg.eigh(a)[1] # eigenvectors
	c = np.diag(1./np.sqrt(lamda + 1)) 
	Sa = np.matrix(S) * (np.matrix(U) * np.matrix(c) * np.matrix(U.T)) # update S every time, scale of ensemble perturbation after the assimilation
	a = np.matrix(yo).T- np.matrix(H) * np.matrix(xmean).T
	K = np.matrix(S) * (np.matrix(U) * np.matrix(np.diag(1./(lamda + 1))) * np.matrix(U.T)) * HS.T * inv(R) # Karman gain for ensemble
	x_mean_assim = np.array(np.matrix(xmean).T + np.matrix(K) * a).flatten()
	x_ens_assim = np.zeros((2,Nens))
	for i in range(len(x_ens_assim)):
		x_ens_assim[:,i] = np.array(np.matrix(xmean).T + np.sqrt(Nens - 1) * Sa[:,i]).flatten()
	return x_ens_assim


	

"""
	

# Compute divergence

for i in range(1,n_it):
	x_pos[:,i] = rgk(x_pos[:,i-1], i, dt) 
	x_noise[:,i] = x_pos[:,i] + np.sqrt(R[1,1]) * np.random.randn(2)
	x_free[:,i] = rgk(x_free[:,i-1], i, dt)
	for k in range(Nens):
		x_ens[:,k,i] = rgk(x_ens[:,k,i-1], i, dt)
		#MP = tanmod_matrix(x_ens[:,k,i-1],P[:,:,i-1],i,dt) 
		#MPT = MP.T
		#P[:,:,i] = tanmod_matrix(x_assim[:,i-1],MPT,i, dt) # P11 - error variance of longitude, P22 - error variance of atitude
		# A POSTERIORI CALCULATIONS 
		#if i % 40 == 0:
		#	K[:,:,i] = np.matrix(P[:,:,i]) * np.matrix(H) * inv(np.matrix(H) * np.matrix(P[:,:,i]) * np.matrix(H.T) + R)
		#	a = np.matrix(x_noise[:,i]).T- np.matrix(H) * np.matrix(x_assim[:,i]).T
		#	x_assim[:,i] = np.array(np.matrix(x_assim[:,i]).T + np.matrix(K[:,:,i]) * a).flatten()
		#	P[:,:,i] = P[:,:,i] - np.matrix(K[:,:,i]) * np.matrix(H) * np.matrix(P[:,:,i])#; break

	
"""	


