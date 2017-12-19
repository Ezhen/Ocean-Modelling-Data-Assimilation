import numpy as np
import matplotlib.pyplot as plt
import sys,os,shutil

num=4000
dx1=1000 #m
h=500. #m
#L=[4000,10000,20000] #m
L=20000. #m if r=2: if L=1000 E=495404968363 dt=431; if L=4000, E= 495404998334; if L=10000, E=495405000323.0 dt=413, if L=20000, E=495405000187 dt=395
#r=[2,3,5] if L=4000 and r=2 dt=424, E= 495404998334; r=3, E= 497857518151.0 dt=424; r=5, E=502762571612 dt=429;
r=2
L1=200000 #m
L2=400000 #m
x=300000 #m
A=1
dt=10 #s
g=9.81
folder = os.path.abspath("Exam")

x=[]; x.append((list(np.arange(0,201000,1000))))
for i in range(1000):
	if x[0][-1]<L2:
		x[0].append(x[0][-1]+r*dx1)
x=x[0]

dx=np.zeros((len(x)-1)); dx=[x[j+1]-x[j] for j in range(len(dx))]


#dx_2=np.zeros((len(dx))); x_2=np.zeros((len(dx)+1)); x_2[0]=x[0]; x_2[1]=x[1];x_2[-1]=x[-1]; dx_2[0]=dx[0]; dx_2[-1]=dx[-1]
#for i in range(1,len(dx)-1):
#	dx_2[i]=dx_2[i-1]*np.cos(i*[])

nu=np.zeros((num,len(dx))); u=np.zeros((num,len(x)))

for j in range(len(dx)):
	nu[0,j]=A*np.exp(-1*(x[j]/L)**2)

#E=sum([(g*nu[0,j]**2+u[0,j]**2/h)*dx[j] for j in range(0,np.where(np.array(x)==L1)[0][0]+1)])
#print E


fig, ax = plt.subplots(figsize=(16, 10))
ax.set_title('Wave propagation',family='Courier New, monospace',fontsize=20, y=0.88)
plt.xlim(0,400000)
plt.ylim(-1,2)
plt.xlabel('Shallow water equation', fontsize=16)
plt.ylabel('Heigth', fontsize=16)

mm=0
for i in range(num-1):				#time loop
	for j in range(len(nu.T)-1):			#space loop
		nu[i+1,j]=nu[i,j]-dt*(u[i,j+1]-u[i,j])/dx[j]
	for j in range(len(u.T)-2):
		u[i+1,j+1]=u[i,j+1]-dt*g*h*(1./dx[j])*(nu[i+1,j+1]-nu[i+1,j])
	print i, sum([(g*nu[i+1,j]**2+u[i+1,j]**2/h)*dx[j] for j in range(0,np.where(np.array(x)==L1)[0][0]+1)])
	#print sum(nu[i+1,:]), sum(u[i+1,:])
	if mm==0 and abs(nu[i+1,np.where(np.array(x)>300000)[0][0]])==nu[i+1,:].max():
		E=sum([(g*nu[i+1,j]**2+u[i+1,j]**2/h)*dx[j] for j in range(0,np.where(np.array(x)==L1)[0][0]+1)])
		mm=mm+1; print E,i+1
	if i%20==0:
		lines=[]
		lines.append(ax.plot([(x[kk+1]-x[kk])/2.+x[kk] for kk in range(len(dx))],nu[i,:], 'b'))
   		ax.lines
		file_name = os.path.abspath(folder+"/tmp"+str(i/20)+".png")
		ax.annotate('time: %s s' %(format((i*20),'.2f')), xy=(0,0),  xycoords='figure fraction',xytext=(910, 100), family='Courier New, monospace',textcoords='offset points',ha="left", va="bottom",fontsize=24, bbox=dict(facecolor='lightgrey', edgecolor='none', pad=5.0))
		fig.savefig(file_name, dpi=100)
   		ax.lines.pop(0)
		print i
		

file_name = os.path.abspath("Exam//tmp%d.png")
os.system("ffmpeg -framerate 12/1 -f image2 -y -i "+file_name+" -r 24 -bit_rate 1800 -vb 20M L_4_r_2.mpeg")
		

