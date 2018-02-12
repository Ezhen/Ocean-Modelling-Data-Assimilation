"""""""""""""""""""""""""""""""""""""""""""""""
"""""" Made by Evgeny Ivanov, 27.01.2018 """"""
"""""""""""""""""""""""""""""""""""""""""""""""

import numpy as np; import matplotlib.pyplot as plt; import sys,os,shutil; from numpy.linalg import inv

def timeoutput(seconds):		# Output time in a good format
	m, s = divmod(seconds, 60)
	h, m = divmod(m, 60)
	return h, m, s

# Function which plots the state-of-art
def plotpng(ii,h,v,ht,vt,hn,vn):			
	lines = []
	lines.append(ax.plot(xx,h, 'b', linewidth = 2, label='assimilated run')); lines.append(ax.plot(xx,ht, 'r', label='true run')); lines.append(ax.plot(xx,hn, 'g', label='free run')); ax.lines
	qui = ax.quiver(xx[0:],h,v[0:-1],np.zeros((num-1)), scale=10, width=0.003, headwidth=1, headlength=1, color='b')
	quit = ax.quiver(xx[0:],ht,vt[0:-1],np.zeros((num-1)), scale=10, width=0.003, headwidth=1, headlength=1, color='r')
	file_name = os.path.abspath(folder+"tmp"+str(ii/out)+".png")
	ax.annotate("time: %d:%02d:%02d" %(timeoutput(ii*dt)), xy=(0,0),  xycoords='figure fraction',xytext=(860, 100), family='Courier New, monospace',textcoords='offset points',ha="left", va="bottom",fontsize=24, bbox=dict(facecolor='lightgrey', edgecolor='none', pad=5.0)); ax.legend(loc=2)
	fig.savefig(file_name, dpi=100)
   	ax.lines.pop(0); ax.lines.pop(0); ax.lines.pop(0); qui.remove(); quit.remove()


folder = "Optimization/"	# the output folder

# Numerical parameters
out = 100			# output and plot each n time step
t = 5000 			# [-] number of time steps
eps = 1e-4
each = 10 #5
##################################################################
##################  Choose the time step  ########################
##################################################################
dt = 0.1	 			# [s] a time step
#dt = 1					# [s] a time step

# Physical constants #
g = 9.81			# [m/s2] the gravity acceleration constant
b = 1000. 			# [m] a wave parameter
a = 2000.			# [m] an a parameter


# Domain parameters  #
x = 10000. 					# [m] the domain's length
num = 101					# [-] descretization of the domain
h = 100. 					# [m] the depth across the domain
dx = x / (num-1)				# [m] the length of one part of the descritized domain

##################################################################
#################  Choose the system noise #######################
##################################################################
R = np.eye(each)/20.				# [cm] uncertainty in the observations = 5 cm
#R = np.eye(each)/10.				# [cm] uncertainty in the observations = 10 cm

# Model initialization
nu_free = np.zeros((t,num-1)) 			# free run WITH assimilation, heights
u_free = np.zeros((t,num))			# free run WITHOUT assimilation, speeds

nu_free_noass = np.zeros((t,num-1)) 			# free run WITHOUT assimilation, heights
u_free_noass = np.zeros((t,num))			# free run WITHOUT assimilation, speeds

nu_true = np.zeros((t,num-1)) 			# true run, heights
nu_noise = np.zeros((t,num-1))			# true run with noise, heights
u_true = np.zeros((t,num))			# true run, speeds

assim = np.zeros((t,num+num-1))			# array containing assimilated heights and speeds

xx = np.zeros((num-1))				# creation of an empty array for distances from the beginning of the domain		

P = np.zeros((num+num-1,num+num-1,t))		# estimation error covarience matrix creation
P_init = np.eye(201); P_init[100:]=0
P[:,:,0] = P_init

K = np.zeros((num+num-1,each,t))
H = np.zeros((each,201))

##################################################################
######### Choose how often you assimilate observations ##########
##################################################################
H[0,0],H[1,10],H[2,20],H[3,30],H[4,40],H[5,50],H[6,60],H[7,70],H[8,80],H[9,90] = 1,1,1,1,1,1,1,1,1,1	# every 10th
#H[0,0],H[1,20],H[2,40],H[3,60],H[4,80] = 1,1,1,1,1							# every 20th

for j in range(num-2):
	xx[j+1] = xx[j] + dx

for j in range(num-1):
	nu_free[0,j] = np.exp(((xx[j]/b)**2)*(-1))			# a loop to initialize the shape
	nu_free_noass[0,j] = np.exp(((xx[j]/b)**2)*(-1))		# of each wave of all three runs
	nu_true[0,j] = np.exp(-1*((xx[j]-a)**2)/(b**2))			# with the initial speed is zero

assim[0] = np.concatenate((nu_free[0,:],u_free[0,:]))


# Create a frame
fig, ax = plt.subplots(figsize=(16, 10))
ax.set_title('Assimilation of shallow water equations',family='Courier New, monospace',fontsize=20, y=0.88)
plt.xlim(0,10000)
plt.ylim(-0.25,1.25)
plt.xlabel('Shallow water equation', fontsize=16)
plt.ylabel('Heigth', fontsize=16)


# descretization scheme (1-order Euler forward) #
def shallow_rev(nu,u):
	nu_new, u_new = np.zeros((100)),np.zeros((101))
	for j in range(num-2):
		u_new[j+1]  = u[j+1]  - dt * g * (nu[j+1]  - nu[j])  /dx
	for j in range(num-1):	
		nu_new[j]  = nu[j]  - dt * h * (u_new[j+1]  - u_new[j])  /dx	
	return nu_new, u_new


# derive tangent-linear model
def tanmod_vector(nu,u,dxx):
	nu_new_1, u_new_1 = shallow_rev(nu+(dxx[0:100]*eps),u)
	nu_new_2, u_new_2 = shallow_rev(nu-(dxx[0:100]*eps),u)
	assim_1 = np.concatenate((nu_new_1,u_new_1))
	assim_2 = np.concatenate((nu_new_2,u_new_2))
	diff = assim_1 - assim_2
	v = diff * 0.5 / eps
	return v

# Than more the flow DIVERGE, than larger the error
def tanmod_matrix(nu, u, P):
	MP = np.zeros((201,201))
	for k in range(201):
		MP[k,:] = tanmod_vector(nu,u,P[k,:])		# product between the tangent linear model and the covariacne matrix P
	return MP

def noise_ass(nnoise):
	noise = np.zeros((each))
	for k in range(each):
		noise[k] = nnoise[k*100/each]
	return noise


# The main function's body (the shallow water equations)
for i in range(1,t):	
	print i		
	nu_free[i],u_free[i] = shallow_rev(nu_free[i-1],u_free[i-1])
	nu_free_noass[i],u_free_noass[i] = shallow_rev(nu_free_noass[i-1],u_free_noass[i-1])
	nu_true[i],u_true[i] = shallow_rev(nu_true[i-1],u_true[i-1])
	nu_noise[i,:] = nu_true[i,:] + R[0,0] * np.random.randn(num-1)#	; print nu_free[i,50]
	assim[i] = np.concatenate((nu_free[i,:],u_free[i,:]))
	noise = noise_ass(nu_noise[i,:])
	MP = tanmod_matrix(nu_free[i-1,:],u_free[i-1,:],P[:,:,i-1])
	MPT = MP.T
	P[:,:,i] = tanmod_matrix(nu_free[i-1,:],u_free[i-1,:],MPT)	# P11 - error variance of longitude, P22 - error variance of atitude
	if i % out == 0 and i>0:
		K[:,:,i] = np.matrix(P[:,:,i]) * np.matrix(H).T * inv(np.matrix(H) * np.matrix(P[:,:,i]) * np.matrix(H).T + R)
		aa = np.matrix(noise).T - (np.matrix(H) * np.matrix(assim[i,:]).T)
		assim[i,:] = (np.matrix(assim[i,:]).T + np.matrix(K[:,:,i]) * np.matrix(aa)).flatten()
		P[:,:,i] = P[:,:,i] - np.matrix(K[:,:,i]) * np.matrix(H) * np.matrix(P[:,:,i])#; break
		nu_free[i],u_free[i] = assim[i,0:100],assim[i,100:]; print i
	if (i-1)%out==0:
		plotpng(i,nu_free[i,:],u_free[i,:],nu_noise[i,:],u_true[i,:],nu_free_noass[i,:],u_free_noass[i,:])

# Copy - paste into an ipython window to plot either the Kalman gain over time, either the covariance matrix over time #
"""
# plot Kalman gain
KK, PP = np.zeros((5000)), np.zeros((5000))
for i in range(len(KK)):
    KK[i]=np.sum(K[:,:,i]); PP[i]=np.sum(P[:,:,i])
fig.clf()
plot(np.arange(t),KK)
ylabel("Kalman gain")
xlabel("Time step")
plot(np.arange(5000),PP)
ylabel("Error covariance matrix")
xlabel("Time step")
"""

# Uncomment if you want to make animation out of assimilation frames, requeires ffmpeg being installed #
#folder = "Opt_good_run_s_01_5cm_each_10/"
#file_name = os.path.abspath(folder+"tmp%d.png")
#os.system("ffmpeg -framerate 2/1 -f image2 -y -i "+file_name+" -r 24 -bit_rate 1800 -vb 20M Anime.mpeg")
