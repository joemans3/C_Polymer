import numpy as np
import sys
import os
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d import axes3d as ax3d
import matplotlib.animation
from cmpaircorrelation import *
from scipy.optimize import curve_fit

import pdb
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


def MSD_tavg(x,y,z,f,N,f_inc = False):
    
    dists = np.zeros(len(x)-1)
    for i in range(len(x)-1):
        dists[i] = dis(x[i],y[i],z[i],x[i+1],y[i+1],z[i+1],N)
    if f_inc == True:
        return np.mean((np.diff(dists/np.diff(f)))**2)/4.
    else:
        return np.mean((np.diff(dists))**2)/6.

def MSD_tavg1(x,y,z,N):
    return np.mean(np.diff(dis(np.array(x)[1:],np.array(y)[1:],np.array(z)[1:],np.array(x)[0],np.array(y)[0],np.array(z)[0],N))**2)/6.


def track_decomp(x,y,z,f,N):
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
        try:
            n_f = np.array(f)[::i]
            msd.append(MSD_tavg(n_x,n_y,n_z,n_f,N))
        except:
            msd.append(MSD_tavg(n_x,n_y,n_z,f,N))

    
    #popt , pcov = curve_fit(fit_MSD,tau,np.array(msd),p0=[1,1],maxfev=10000)
    
    
    return np.array(msd)




def track_decomp1(x,y,f):
    #takes tracks and finds MSD for various timestep conditions.
    
    #return array-like:
    #msd = msd values at all tau values considered
    #popt = fitted parameters on MSD equation
    #pcov = covariance matrix of fit
    max_track_decomp = 1.
    max_decomp = np.floor(len(x)/max_track_decomp)
    tau = range(1,int(max_decomp+1.0))
    msd = []
    for i in tau:
        if i < len(x):
            n_x = np.array(x)[::i]
            n_y = np.array(y)[::i]
            try:
                n_f = np.array(f)[::i]
                msd.append(MSD_tavg(n_x,n_y,n_f))
            except:
                msd.append(MSD_tavg(n_x,n_y,f))

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
            VA_data_type.append(np.resize(data[T_P*j + sizeMP[i-1]:sizeMP[i]+T_P*j],(sizeP[i],sizeM[i],4)))
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


def fit_decom(d,thresh = 0.5,max = 100000):

    popt, pcov = curve_fit(fit_MSD,range(1,len(d)+1)[:int(len(range(1,len(d)+1))*thresh)],d[:int(len(range(1,len(d)+1))*thresh)],p0=[1,1],maxfev=10000)
    return [popt,pcov]




############################################################################################
#msd stuff for ALL POLYMER SIMULATIONS

new_size = np.asarray(sizeM)*np.asarray(sizeP)
#tracks per polymer type
track_ptype = []


msd_t_meana = np.zeros(len(sizeM))
msd_t_mean_verba = []
    
fit_verbosea = [[] for i in new_size]
fit_verbosed = [[] for i in new_size]
fit_msdaa = []
    
msd_s_meana = []
msd_s_mean_va = [[] for i in new_size]

for k in range(len(sizeM)):
    n_ordered_arr = np.zeros((sizeM[k],sizeP[k],samples[0],4))
    for i in range(len(VA_data_type_O)):
        for j in range(len(VA_data_type_O[i][k])):
            for l in range(len(VA_data_type_O[i][k][j])):
                n_ordered_arr[l][j][i] = VA_data_type_O[i][k][j][l]

    track_ptype.append(n_ordered_arr)

    #msd

    msd_t_mean = np.zeros(len(n_ordered_arr))


    fit_verbose = []

    msd_s_mean = []
    msd_s_mean_v = []
    for i in range(len(n_ordered_arr)):
        msd_t_mean_p = np.zeros(len(n_ordered_arr[i]))
        msd_dec = []
        
        for j in range(len(n_ordered_arr[i])):
            pp = n_ordered_arr[i][j]
            
            msd_t_mean_p[j] = MSD_tavg(pp[:,0],pp[:,1],pp[:,2],1,sizeN)
            tttt=track_decomp(pp[:,0],pp[:,1],pp[:,2],1,sizeN)
            msd_dec.append(tttt)
            msd_s_mean_va[k].append(tttt)
            aaaa=fit_decom(tttt)
            fit_verbosea[k] = fit_verbosea[k] + [aaaa[0][1]]
            fit_verbosed[k] = fit_verbosed[k] + [aaaa[0][0]]
        


        msd_dec = np.mean(np.asarray(msd_dec),axis = 0)
        msd_s_mean.append(msd_dec)
                        
        msd_t_mean[i] = np.mean(msd_t_mean_p)

    msd_s_meana.append(np.mean(np.asarray(msd_s_mean),axis = 0))

    fit_msdaa.append(fit_decom(np.mean(np.asarray(msd_s_mean),axis = 0))[0])

    msd_t_meana[k] = np.mean(msd_t_mean)
    msd_t_mean_verba.append(msd_t_mean)




def msd_plot(m,verbose = False):
    '''Verbose: all msd or averaged ; True, False'''
    if verbose:
        for i in range(len(m)):
            for j in m[i]:
                plt.plot(j)
            plt.ylabel("MSD")
            plt.xlabel("Tau")
            plt.title("P = {0}, M = {1}".format(sizeP[i],sizeM[i]))
            plt.yscale("log")
            plt.xscale("log")
            plt.show()
    else:
        for i in range(len(m)):
            plt.plot(m[i],label="{0}".format(i+1))
        plt.ylabel("MSD")
        plt.xlabel("Tau")
        plt.title("P = {0}, M = {1}".format(sizeP,sizeM))
        plt.yscale("log")
        plt.xscale("log")
        plt.legend()
        plt.show()
    return


def create_box_plot(box_data,tick_list,y_label = "",x_label = "",y_lim = (),title = ""):
    ticks = tick_list
    plt.boxplot(box_data,positions = range(1,len(tick_list)+1))
    for i in range(1,len(tick_list)+1):
        y = box_data[i-1]
        x = np.random.normal(i, 0.04, size=len(y))
        plt.plot(x, y, 'r.', alpha=0.2)
    try:
        plt.ylim(y_lim)
    except:
        print "Warning: y_lim not valid"
    plt.xticks(xrange(1, len(ticks) * 1 + 1, 1), ticks)
    plt.ylabel(y_label)
    plt.xlabel(x_label)
    plt.title(title)
    plt.show()

    return

from sklearn import mixture

def GMM_utility(data, n, biners=50, inclusion_thresh = [0,100], verbose=True, title_1d="", title_2d="", x_label="", y_label_2d="", log=True, x_limit = ()):
    
    data = np.array(data)
    
    p_thresh = np.percentile(data,inclusion_thresh)
    inds = ((data>=p_thresh[0]) & (data<=p_thresh[1]))
    data = data[inds]
    
    gmix = mixture.GMM(n_components=n, covariance_type='diag')
    if log:
        (results,bins) = np.histogram(np.log(data),density='true',bins=biners)
    else:
        (results,bins) = np.histogram(data,density='true',bins=biners)
    
    
    data_arr = np.zeros((len(data),2))
    data_arr[:,0] = np.random.normal(1, 0.04, size=len(data))
    if log:
        data_arr[:,1] = np.log(data)
    else:
        data_arr[:,1] = data
    if verbose:
        plt.plot(data_arr[:,1],data_arr[:,0],'r.')
        plt.ylim((0,2))
        plt.title(title_1d)
        plt.xlabel(x_label)
        plt.show()
    gmix.fit(data_arr)

    if log:
        print "Fitted Mean: {0} +/- {1}".format(gmix.means_[:,1],np.sqrt(gmix.covars_[:,1]))
        print "Fitted Mean(normal): {0} +/- {1}".format(np.exp(gmix.means_[:,1]),np.exp(gmix.means_[:,1])*np.sqrt(gmix.covars_[:,1]))
    else:
        print "Fitted Mean: {0} +/- {1}".format(gmix.means_[:,1],np.sqrt(gmix.covars_[:,1]))
    max_r = np.max(results)
    plt.plot(np.diff(bins)+bins[:len(bins)-1],results)

    for i in gmix.means_:
        plt.axvline(x=i[1],color='red')
    plt.title(title_2d)
    plt.xlabel(x_label)
    plt.ylabel(y_label_2d)
    try:
        plt.xlim(x_limit)
    except:
        print "Warning: x_limit is invalid"
    plt.show()

    return

























##########################################################################################################
#testing out new animations for line plot only!
##########################################################################################################

#fig = plt.figure()
#ax1 = ax3d.Axes3D(fig)
#line, = ax1.plot([], [], lw=1, linestyle= "-")
#
#ax1.set_xlim(0,sizeN)
#ax1.set_ylim(0,sizeN)
#ax1.set_zlim(0,sizeN)
#ax1.set_xlabel("x")
#ax1.set_ylabel("y")
#ax1.set_zlabel("z")
#ax1.text2D(0, 0, "Title", transform=ax1.transAxes)
#
#plotlays, plotcols = [np.sum(sizeP)], ["orange","black"]
#lines = []
#for index in range(np.sum(sizeP)):
#    lobj = ax1.plot([],[],lw=1)[0]
#    lines.append(lobj)
#
#def init():
#    for line in lines:
#        line.set_data([],[])
#        line.set_3d_properties([])
#    return lines
#
#def animate(i):
#
#    
#    xlist = []
#    ylist = []
#    zlist = []
#    
#    for j in range(len(VA_data_type_O[i])):
#        for k in range(len(VA_data_type_O[i][j])):
#            xlist.append(np.asarray(VA_data_type_O[i][j][k])[:,0])
#            ylist.append(np.asarray(VA_data_type_O[i][j][k])[:,1])
#            zlist.append(np.asarray(VA_data_type_O[i][j][k])[:,2])
#    for lnum,line in enumerate(lines):
#        line.set_data(xlist[lnum], ylist[lnum])
#        line.set_3d_properties(zlist[lnum])
#    return lines
#
#anim = matplotlib.animation.FuncAnimation(fig, animate,init_func=init, frames=range(1,samples[0]), blit=True)
#
#plt.show()

##########################################################################################################
##########################################################################################################




































#############################################################################################
#
##plotting and animations
#
#############################################################################################
#
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








