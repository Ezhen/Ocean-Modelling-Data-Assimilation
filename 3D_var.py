import numpy as np
import matplotlib.pyplot as plt

dx = 0.05
x = np.arange(0,2,dx)
f = np.sin(2*np.pi*x)
R=1 # error of our observation
yo = 1
L = 0.5 # a length scale of the cost function


H = np.zeros((1,40)); H=np.matrix(H); H[0,9] = 1

def lap(f):
	"""Laplacian operator"""
	laplap = (np.roll(f,1)-2*f+np.roll(f,-1))/(dx**2)
	return laplap

def K(f):
	"""Cost function"""
	KK = L**4 * np.dot(lap(f),lap(f).T) + np.dot(f,f) + (H*np.matrix(f).T - yo).T * R * (H*np.matrix(f).T - yo)
	return KK

def grad(f):
	s = H*np.matrix(f).T - yo
	grd = 2 * (L**4 * lap(lap(f)) + f + np.array(H.T * R * s)[:,0])
	return grd

# H*np.matrix(f).T

# df=np.random.randn(len(x))/1000000.
#K(f+df) - K(f)
#np.dot(df.T,grad(f))

fa = np.zeros((len(f)))
Nint=200000
a=0.0000001
Ji = np.zeros((Nint))
for i in range(Nint):
	fa += - a * grad(fa)
	if i%1000==0:
		print i,K(fa)
	#Ji[i] = K(fa) # outputiing a cost function to keep track of the convergence
plt.plot(x,fa)
plt.show()
"""

der1 = (f[1:]-f[0:-1])/dx
#der2 = (der1[1:]-der1[0:-1])/dx
#der2 = (f[2:]-2*f[1:-1]+f[0:-2])/(dx**2)
der2 = lap(f)

print np.shape(f),np.shape(der1),np.shape(der2)
plt.plot(x,f,label='f')
#plt.plot(x[0:-1]+dx/2.,der1,label='der1')
plt.plot(x,der2,label='der2')
plt.legend(loc=2)
plt.show()
"""
"""
# Laplacian test
a=np.random.randn(100)
b=np.random.randn(100)
c1 = np.dot(a.T,lap(b))
c2 = np.dot(b.T,lap(a))
print c1,c2
"""
