import numpy as np
import sys
import os
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation
from cmpaircorrelation import *
path = raw_input("Path of folder: ")





####calculating center of mass with periodic conditions

def cm(x,y,z,sizeN):
    #transform x,y to -pi <-> pi
    xpi=x*2.*np.pi/sizeN
    ypi=y*2.*np.pi/sizeN
    zpi=z*2.*np.pi/sizeN
    #find the geometric mean (all points have weighting factor of 1)
    xpi_meanc=np.mean(np.cos(xpi))
    xpi_means=np.mean(np.sin(xpi))
    
    ypi_meanc=np.mean(np.cos(ypi))
    ypi_means=np.mean(np.sin(ypi))
    
    
    zpi_meanc=np.mean(np.cos(zpi))
    zpi_means=np.mean(np.sin(zpi))
    
    
    #transform back to x,y space
    thetax=np.arctan2(-xpi_means,-xpi_meanc) + np.pi
        
    thetay=np.arctan2(-ypi_means,-ypi_meanc) + np.pi
    
    thetaz=np.arctan2(-zpi_means,-zpi_meanc) + np.pi
        
    xcm=sizeN*thetax/(2.*np.pi)
    ycm=sizeN*thetay/(2.*np.pi)
    zcm=sizeN*thetaz/(2.*np.pi)
    
    return np.array([xcm,ycm,zcm])







#system variables
sizeM = []
sizeP = []
sizeN = 0
samples = 0

def convert_to_int(input):

    return map(int,input.split(","))

def convert_globals(input):
    SC_index_a = []
    C_index_a = []
    sizeN = 0
    sizeM = []
    sizeP = []
    for i in range(0,len(input)):
        if input[i]==";":
            SC_index_a.append(i)
        if input[i]==",":
            C_index_a.append(i)
    sizeN = int(input[:SC_index_a[0]])
    sizeP = convert_to_int(input[SC_index_a[0]+1:SC_index_a[1]])
    sizeM = convert_to_int(input[SC_index_a[1]+1:SC_index_a[2]])
    samples = convert_to_int(input[SC_index_a[2]+1:])
    return [sizeN, sizeP, sizeM, samples]


with open(path,"r") as f:
    line = (f.readline()).rstrip('\n')
    sizeN,sizeP,sizeM,samples = convert_globals(line)
f.close()

sizeMP = np.cumsum(np.asarray(sizeM)*np.asarray(sizeP))
T_P= np.sum(np.asarray(sizeM)*np.asarray(sizeP))
data = np.loadtxt(path,delimiter=",",skiprows=1)


VA_data_type_O = []  #datatype for line plots, doesnt work for animation

for j in range(0,samples[0]):

    VA_data_type = []

    for i in range(0,len(sizeMP)):
        if i==0:
            VA_data_type.append(np.resize(data[T_P*j:sizeMP[i]+T_P*j],(sizeP[i],sizeM[i],4)))
        else:
            VA_data_type.append(np.resize(data[T_P*j + sizeMP[i-1]+1:sizeMP[i]+T_P*j],(sizeP[i],sizeM[i],4)))
    VA_data_type_O.append(VA_data_type)

#finding the center of mass for each frame

CM_ar=[]
for i in VA_data_type_O:
    CM_ar.append(cm(i[0][:,0],i[0][:,1],i[0][:,2],sizeN))




#calculating pair correlation function.

pc_holder=[]
radius_holder=[]
radius_holder2=[]
for i in range(len(VA_data_type_O)):
    temp1 , temp2 = paircorrelation3D(VA_data_type_O[i][0][:,0],VA_data_type_O[i][0][:,1],VA_data_type_O[i][0][:,2],sizeN,CM_ar[i],VA_data_type_O[i][0][:,3],dr=0.5)
    pc_holder.append(temp1)
    radius_holder.append(temp2)
#radius_holder2.append(temp3)
plt.plot(1./(np.array(radius_holder)[:,1]))
plt.title("Correlation Lenght Fit With Exponential Only")
plt.xlabel("Sample Frame")
plt.ylabel("Correlation Lenght (units of lattice)")
plt.show()

for i in pc_holder:
    plt.plot(i,'ro')
    plt.yscale("log")
    plt.show()
'''
plt.plot(1./(np.array(radius_holder2)[:,1]))
plt.title("Correlation Lenght Fit With Exponential + Power Law Fit")
plt.xlabel("Sample Frame")
plt.ylabel("Correlation Lenght (units of lattice)")
plt.show()
'''

#Creating datatype for animation (does not include animation)
'''
animate_DT = []

for j in range(0,samples[0]):
    animate_DT.append(np.reshape(VA_data_type_O[j],(T_P,4)))
'''


temp=[]

for j in VA_data_type_O:
    tmp=[]
    for i in j:
        for k in i:
            for uu in k:
                tmp.append(uu)
    temp.append(np.array(tmp))

animate_DT=temp





############################################################################################

#plotting and animations

############################################################################################

# Set up formatting for the movie files
Writer = matplotlib.animation.writers['ffmpeg']
writer = Writer(fps=10, metadata=dict(artist='Baljyot Parmar, all rights reserved'), bitrate=1800)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
title = ax.set_title('3D Test')
def update_graph(num):
    i = animate_DT[num]

    graph._offsets3d = (i[:,0],i[:,1],i[:,2])

    title.set_text('3D Test, Sample={}'.format(num))

graph = ax.scatter(animate_DT[0][:,0],animate_DT[0][:,1],animate_DT[0][:,2],c=animate_DT[0][:,3])

ax.set_xbound(lower=-1,upper=sizeN+1)
ax.set_ybound(lower=-1,upper=sizeN+1)
ax.set_zbound(lower=-1,upper=sizeN+1)
ani = matplotlib.animation.FuncAnimation(fig, update_graph, range(1,samples[0]), blit=False)
##ani.save('/Users/baljyot/Documents/Polymer_Output/new_ani.mp4',writer=writer)
ani.save('{0}/new_ani.mp4'.format(os.getcwd()),writer=writer)
plt.show()


'''
for k in range(0,samples[0]):
    for i in VA_data_type_O[k]:
        for j in i:
            ax.scatter(j[:,0],j[:,1],j[:,2],c=j[:,3])
    #ax.plot(j[:,0],j[:,1],j[:,2],'-b')
    ax.set_xbound(lower=0,upper=sizeN)
    ax.set_ybound(lower=0,upper=sizeN)
    ax.set_zbound(lower=0,upper=sizeN)
plt.show()


'''








