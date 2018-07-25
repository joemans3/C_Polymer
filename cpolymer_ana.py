import numpy as np
import sys
import os
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation
from cmpaircorrelation import *
from scipy.optimize import curve_fit
path = raw_input("Path of folder: ")




def dis(x,y,z,x1,y1,z1,N):
    t1=np.minimum(abs(x1-x),abs(x1-x - N))
    tx=np.minimum(t1,abs(x1-x + N))
    t2=np.minimum(abs(y1-y),abs(y1-y - N))
    ty=np.minimum(t2,abs(y1-y + N))
    t3=np.minimum(abs(z1-z),abs(z1-z - N))
    tz=np.minimum(t3,abs(z1-z + N))
    temp=np.sqrt((tx)**2 + (ty)**2 + (tz)**2)
    return temp


def fit_MSD(t,p_0,p_1):
    return p_0 * (t**(p_1))


def MSD_tavg(x,y,z,N):
    return np.mean(np.diff(dis(np.array(x)[1:],np.array(y)[1:],np.array(z)[1:],np.array(x)[0],np.array(y)[0],np.array(z)[0],N))**2)/6.


def track_decomp(x,y,z,N):
    #takes tracks and finds MSD for various timestep conditions.
    
    #return array-like:
    #msd = msd values at all tau values considered
    #popt = fitted parameters on MSD equation
    #pcov = covariance matrix of fit
    max_track_decomp = 10.
    max_decomp = np.floor(len(x)/max_track_decomp)
    tau = range(1,int(max_decomp+1))
    msd = []
    for i in tau:
        n_x = np.array(x)[::i]
        n_y = np.array(y)[::i]
        n_z = np.array(z)[::i]
        msd.append(MSD_tavg(n_x,n_y,n_z,N))
    
    #popt , pcov = curve_fit(fit_MSD,tau,np.array(msd),p0=[1,1],maxfev=10000)
    
    
    return np.array(msd)























def epo1(x,p1,p2):
    return p1*np.exp(-x*p2)
def epo(x,*p0):
    return (1./(x**p0[0]))*p0[1]*np.exp(-x*p0[2]) #+ p0[3]

def dist(x,y,z,c,N):
    t1=np.minimum(abs(c[0]-x),abs(c[0]-x - N))
    tx=np.minimum(t1,abs(c[0]-x + N))
    t2=np.minimum(abs(c[1]-y),abs(c[1]-y - N))
    ty=np.minimum(t2,abs(c[1]-y + N))
    t3=np.minimum(abs(c[2]-z),abs(c[2]-z - N))
    tz=np.minimum(t3,abs(c[2]-z + N))
    temp=np.sqrt((tx)**2 + (ty)**2 + (tz)**2)
    return temp


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
#
#CM_ar=[]
#x_per_frame=[]
#y_per_frame=[]
#z_per_frame=[]
#c_per_frame=[]
#for i in VA_data_type_O:
#    per_frame_per_p_x =[]
#    per_frame_per_p_y =[]
#    per_frame_per_p_z =[]
#    per_frame_per_p_c =[]
#    for j in i[0]:
#        per_frame_per_p_x.append(j[:,0])
#        per_frame_per_p_y.append(j[:,1])
#        per_frame_per_p_z.append(j[:,2])
#        per_frame_per_p_c.append(j[:,3])
#    fx=np.array(per_frame_per_p_x).flatten()
#    fy=np.array(per_frame_per_p_y).flatten()
#    fz=np.array(per_frame_per_p_z).flatten()
#    fc=np.array(per_frame_per_p_c).flatten()
#    x_per_frame.append(fx)
#    y_per_frame.append(fy)
#    z_per_frame.append(fz)
#    c_per_frame.append(fc)
#    CM_ar.append(cm(fx,fy,fz,sizeN))
#
#
###########
##finding total distance per frame
#distance_arr=[]
#for i in range(len(x_per_frame)):
#    distance_t=0
#    for j in range(len(x_per_frame[i])):
#        for kk in range(j+1,len(x_per_frame[i])):
#            distance_t+=dist(x_per_frame[i][j],y_per_frame[i][j],z_per_frame[i][j],[x_per_frame[i][kk],y_per_frame[i][kk],z_per_frame[i][kk]],sizeN)
#
#    distance_arr.append(distance_t)
#
#plt.plot(distance_arr)
#plt.title("Total Distance per Sample")
#plt.xlabel("Sample #")
#plt.ylabel("Distance (lattice units)")
#plt.show()





'''

#calculating pair correlation function.

pc_holder=[]
radius_holder=[]
radius_holder2=[]
distance=[]
radi=[]
popt1=[]

for i in range(len(VA_data_type_O)):
    temp1 , temp2 ,temp3 , temp4 ,temp5 = paircorrelation3D(x_per_frame[i],y_per_frame[i],z_per_frame[i],sizeN,CM_ar[i],c_per_frame[i],dr=2)
    pc_holder.append(temp1)
    radius_holder.append(temp2)
    distance.append(temp3)
    radi.append(temp4)
    popt1.append(temp5)

#radius_holder2.append(temp3)
plt.plot(1./(np.array(radius_holder)[:,1]))
plt.title("Correlation Lenght Fit With Exponential Only")
plt.xlabel("Sample Frame")
plt.ylabel("Correlation Lenght (units of lattice)")
plt.show()


plt.plot((np.array(popt1)[:,0]))
plt.title("Correlation Lenght Fit With Exponential Only")
plt.xlabel("Sample Frame")
plt.ylabel("Correlation Lenght (units of lattice)")
plt.show()

'''




##############################################################################################
####correlation length stuff


#
#pc_holder=[]
#radius_holder=[]
#radius_holder2=[]
#radi=[]
#popt1=[]
#
#for i in range(len(VA_data_type_O)):
#    temp1, temp2, temp3 = paircorrelation3D_a(x_per_frame[i],y_per_frame[i],z_per_frame[i],sizeN,CM_ar[i],c_per_frame[i],dr=0.5)
#    pc_holder.append(temp1)
#    radius_holder.append(temp2)
#    radi.append(temp3)
#
#'''
##radius_holder2.append(temp3)
#for i in range(len(pc_holder)):
#    plt.plot(radi[i],pc_holder[i],'ro')
#    plt.plot(radi[i],epo1(radi[i],radius_holder[i][0],radius_holder[i][1]))
#    plt.title("Correlation Function")
#    plt.xlabel("units")
#    plt.ylabel("Averaged Correlation Function")
#    plt.yscale("log")
#    #plt.xscale("log")
#    plt.show()
#
#'''
#
#plt.plot(1./(np.array(radius_holder)[:,1]))
#plt.title("Correlation Length w Exponential Over Samples")
#plt.xlabel("Sample #")
#plt.ylabel("Correlation Length (lattice units)")
#plt.show()
#







#for i in range(len(distance)):
#    plt.hist(pc_holder[i])
#
#    plt.show()


'''
for i in range(len(pc_holder)):
    plt.plot(np.arange(0,sizeN+0.5,0.5),epo(np.arange(0,sizeN+0.5,0.5),popt1[i][0],popt1[i][1],popt1[i][2]))
    plt.plot(np.arange(0,sizeN+0.5,0.5),pc_holder[i],'ro')
    #plt.yscale("log")
    plt.show()
'''
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
#msd stuff for SINGLE POLYMER SIMULATIONS ONLY!




'''

n_ordered_arr = np.zeros((sizeM[0],samples[0],4))
for i in range(len(animate_DT)):
    for j in range(len(animate_DT[i])):
        n_ordered_arr[j][i] = animate_DT[i][j]

mt_msd = np.zeros(sizeM[0])
mt_fit = np.zeros((sizeM[0],2))
for i in range(len(n_ordered_arr)):
    pp = n_ordered_arr[i]
    one_msd_decomp = track_decomp(pp[:,0],pp[:,1],pp[:,2],sizeN)
    popt, pcov = curve_fit(fit_MSD,range(1,len(one_msd_decomp)+1)[:int(len(range(1,len(one_msd_decomp)+1))*0.5)],one_msd_decomp[:int(len(range(1,len(one_msd_decomp)+1))*0.5)],p0=[1,1],maxfev=10000)
    mt_fit[i] = popt
    mt_msd[i] = MSD_tavg(pp[:,0],pp[:,1],pp[:,2],sizeN)

plt.plot(mt_msd)
plt.show()

plt.plot(mt_fit[:,1])
plt.show()

'''




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
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_zlabel("z")
ani = matplotlib.animation.FuncAnimation(fig, update_graph, range(1,samples[0]), blit=False)
##ani.save('/Users/baljyot/Documents/Polymer_Output/new_ani.mp4',writer=writer)
ani.save('{0}/new_ani.mp4'.format(os.getcwd()),writer=writer)
plt.show()

plt.plot(animate_DT[0][:,0],animate_DT[0][:,1],'ro')
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








