import numpy as np
import sys
import os
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation


path = raw_input("Path of folder: ")


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

#Creating datatype for animation (does not include animation)

animate_DT = []

for j in range(0,samples[0]):
    animate_DT.append(np.reshape(VA_data_type_O[j],(T_P,4)))

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

ax.set_xbound(lower=0,upper=sizeN+1)
ax.set_ybound(lower=0,upper=sizeN+1)
ax.set_zbound(lower=0,upper=sizeN+1)
ani = matplotlib.animation.FuncAnimation(fig, update_graph, range(1,samples[0]), blit=False)
ani.save('/Users/baljyot/Documents/Polymer_Output/new_ani.mp4',writer=writer)
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








