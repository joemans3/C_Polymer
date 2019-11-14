import numpy as np
import sys
import os
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d import axes3d as ax3d
import matplotlib.animation
import matplotlib.colors as colors
import math
from pathlib import Path

from cmpaircorrelation import *
from scipy.optimize import curve_fit
import shutil
import pdb

#convert degrees to rad
def d_to_rad(deg_):
    return np.array(deg_)*np.pi/180.0

#convert rad to deg
def rad_to_d(rad_):
    return np.array(rad_)*180.0/np.pi




def dot(a,b):
    return a[0]*b[0] + a[1]*b[1]


def ang(a,b):

    ''' takes input as tuple of tuples of X,Y.'''

    la = [(a[0][0]-a[1][0]), (a[0][1]-a[1][1])]
    lb = [(b[0][0]-b[1][0]), (b[0][1]-b[1][1])]

    dot_ab = dot(la,lb)

    ma = dot(la,la)**0.5
    mb = dot(lb,lb)**0.5

    a_cos = dot_ab/(ma*mb)
    try:
        angle = math.acos(dot_ab/(mb*ma))
    except:
        angle = math.acos(round(dot_ab/(mb*ma)))

    ang_deg = math.degrees(angle)%360

    if ang_deg-180>=0:
        return 360 - ang_deg
    else:
        return ang_deg


def angle_trajectory_3d(x,y,z,ref = True):
    ''' Takes input (x,y,z) of a series of arrays or one array of 
    trajectorie(s) and returns a series or one array of angles in 3D.
    
    INPUTS:

    x,y (array-like): series of arrays or one array of trajectorie(s).

    ref (boolian): If True, return an extra angle which is the angle between the first line in the set and a verticle line
    
    RETURN:

    Array-like: A series or one array of angles in 2D depending on shape of x,y,z.


    '''
   
    if isinstance(x[0],list):
        angle_list = [[] for i in x]
        for i in range(len(x)):
            for j in range(len(x[i])-2):

                angle_list[i].append(ang((((x[i],y[i],z[i]),(x[i+1],y[i+1],z[i+1])),((x[i+1],y[i+1],z[i+1]),(x[i+2],y[i+2],z[i+2])))))
        return angle_list

    else:
        return [ang(((x[i],y[i],z[i]),(x[i+1],y[i+1],z[i+1])),((x[i+1],y[i+1],z[i+1]),(x[i+2],y[i+2],z[i+2]))) for i in range(len(x)-2)]


def dis(x,y,z,x1,y1,z1,N):
    t1=np.minimum(abs(x1-x),abs(x1-x - N))
    tx=np.minimum(t1,abs(x1-x + N))
    t2=np.minimum(abs(y1-y),abs(y1-y - N))
    ty=np.minimum(t2,abs(y1-y + N))
    t3=np.minimum(abs(z1-z),abs(z1-z - N))
    tz=np.minimum(t3,abs(z1-z + N))
    temp=np.sqrt((tx)**2 + (ty)**2 + (tz)**2)
    return temp




#plot a histogram of angles in polar coordinates for each fraction
def hist_polar(data, fig = False, ax_n = False, bin_n = 10, show = True, include_ = True, align = 'edge'): #align can also be 'center'
    '''Helper function to plot histogram of angles in deg on polar coordinates.'''

    bins = np.linspace(0.0, 2.0 * np.pi, bin_n + 1.0)

    #convert to radians
    temp = d_to_rad(data)
    if include_:
        angle_rad = temp
    else:
        angle_rad = temp[temp<3.14]

    print(angle_rad)
    n, _ = np.histogram(angle_rad, bins)

    width = 2.0 * np.pi / bin_n

    if fig == False:
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='polar')

    else:
        ax = fig.add_subplot(ax_n)


    bars = ax.bar(bins[:bin_n], n, width=width, bottom=0.0,align = align)

    for bar in bars:
        bar.set_alpha(0.5)

    if show:
        plt.show()

    return angle_rad




def fit_MSD(t,p_0,p_1):
    return p_0 * (t**(p_1))


def MSD_tavg(x,y,z,f,N,f_inc = False):
    
    dists = np.zeros(len(x)-1)
    for i in range(len(x)-1):
        dists[i] = dis(x[i],y[i],z[i],x[i+1],y[i+1],z[i+1],N)
    if f_inc == True:
        return np.mean((np.diff(dists/np.diff(f)))**2)#/4.
    else:
        return np.mean((np.diff(dists))**2)#/6.

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
    tau = list(range(1,int(max_decomp+1)))
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
    tau = list(range(1,int(max_decomp+1.0)))
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










path = input("Give me path ")







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
R = 0
temp = 0

def convert_to_int(input):

    return list(map(int,input.split(",")))
def convert_to_float(input):
    return list(map(float,input.split(",")))
def convert_globals(input, input2):
    SC_index_a = []
    C_index_a = []
    SC_index2_a = []
    C_index2_a = []
    sizeN = 0
    sizeM = []
    sizeP = []

    for i in range(0,len(input)):
        if input[i]==";":
            SC_index_a.append(i)
        if input[i]==",":
            C_index_a.append(i)

    for i in range(0,len(input2)):
        if input2[i]==",":
            C_index2_a.append(i)
        if input2[i]==":":
            SC_index2_a.append(i)

    R = int(input2[SC_index2_a[0]+1:C_index2_a[0]])
    temp = float(input2[SC_index2_a[2]+1:])
    sizeN = int(input[:SC_index_a[0]])
    sizeP = convert_to_int(input[SC_index_a[0]+1:SC_index_a[1]])
    sizeM = convert_to_int(input[SC_index_a[1]+1:SC_index_a[2]])
    samples = convert_to_int(input[SC_index_a[2]+1:SC_index_a[3]])
    version = input[SC_index_a[3]+1:]
    return [sizeN, sizeP, sizeM, samples,version,R,temp]

#change to ignore first line
with open(path,"r") as f:
    line = 0
    line2 = 0
    for i in range(3):
        temp_line = (f.readline()).rstrip('\n')
        if i == 0:
            line = temp_line
        elif i == 2:
            line2 = temp_line
    sizeN,sizeP,sizeM,samples,version,R,temp = convert_globals(line, line2)
f.close()

dname = os.path.dirname(path)
os.chdir(dname)




if (not os.path.exists(path[len(dname)+1:]+"_folder")) and (str(Path.cwd())[-6:] != "folder"):
    os.makedirs(path[len(dname)+1:]+"_folder")

    os.chdir(path+"_folder")

    shutil.move(path,path+"_folder/"+path[len(dname)+1:])
    os.chdir(os.path.dirname(path+"_folder/"+path[len(dname)+1:]))
    path = path+"_folder/"+path[len(dname)+1:]







sizeMP = np.cumsum(np.asarray(sizeM)*np.asarray(sizeP))
T_P= np.sum(np.asarray(sizeM)*np.asarray(sizeP))
data = np.loadtxt(path,delimiter=",",skiprows=3)

write_bool = False







VA_data_type_O = []  #datatype for line plots, doesnt work for animation
with open(path[len(dname)+1:]+".xyz","w+") as f:

    for j in range(0,samples[0]):
        f.write("{0}\n".format(T_P))
        f.write("\n")
        VA_data_type = []
        if write_bool:
            f.write("H {0} {1} {2}\n".format(data[j][0],data[j][1],data[j][2]))
        for i in range(0,len(sizeMP)):
            if i==0:
                VA_data_type.append(np.resize(data[T_P*j:sizeMP[i]+T_P*j],(sizeP[i],sizeM[i],4)))
            else:
                VA_data_type.append(np.resize(data[T_P*j + sizeMP[i-1]:sizeMP[i]+T_P*j],(sizeP[i],sizeM[i],4)))
        VA_data_type_O.append(VA_data_type)
f.close()
#finding the center of mass for each frame

CM_ar=[]
x_per_frame=[]
y_per_frame=[]
z_per_frame=[]
c_per_frame=[]
for i in VA_data_type_O:
   per_frame_per_p_x =[]
   per_frame_per_p_y =[]
   per_frame_per_p_z =[]
   per_frame_per_p_c =[]
   for j in i[0]:
       per_frame_per_p_x.append(j[:,0])
       per_frame_per_p_y.append(j[:,1])
       per_frame_per_p_z.append(j[:,2])
       per_frame_per_p_c.append(j[:,3])
   fx=np.array(per_frame_per_p_x).flatten()
   fy=np.array(per_frame_per_p_y).flatten()
   fz=np.array(per_frame_per_p_z).flatten()
   fc=np.array(per_frame_per_p_c).flatten()
   x_per_frame.append(fx)
   y_per_frame.append(fy)
   z_per_frame.append(fz)
   c_per_frame.append(fc)
   CM_ar.append(cm(fx,fy,fz,sizeN))





#finding total distance per frame
# distance_arr=[]
# for i in range(len(x_per_frame)):
#    distance_t=0
#    for j in range(len(x_per_frame[i])):
#        for kk in range(j+1,len(x_per_frame[i])):
#            distance_t+=dist(x_per_frame[i][j],y_per_frame[i][j],z_per_frame[i][j],[x_per_frame[i][kk],y_per_frame[i][kk],z_per_frame[i][kk]],sizeN)

#    distance_arr.append(distance_t)

# plt.plot(distance_arr)
# plt.title("Total Distance per Sample")
# plt.xlabel("Sample #")
# plt.ylabel("Distance (lattice units)")
# plt.show()







#calculating pair correlation function.

# pc_holder=[]
# radius_holder=[]
# radius_holder2=[]
# distance=[]
# radi=[]


# for i in range(len(VA_data_type_O)):
#     temp1 , temp2 ,temp3 , temp4 = paircorrelation3D(x_per_frame[i],y_per_frame[i],z_per_frame[i],sizeN,CM_ar[i],c_per_frame[i],dr=1)
#     pc_holder.append(temp1)
#     radius_holder.append(temp2)
#     distance.append(temp3)
#     radi.append(temp4)


# # for i in range(len(pc_holder[::10])):
# #     plt.plot(radi[i],pc_holder[i],'ro')
# #     plt.plot(radi[i],epo1(radi[i],radius_holder[i][0],radius_holder[i][1]))
# #     plt.title("Correlation Function")
# #     plt.xlabel("units")
# #     plt.ylabel("Averaged Correlation Function")
# #     plt.yscale("log")
# #     #plt.xscale("log")
# #     plt.show()

# radius_holder2.append(temp3)
# plt.plot(1./(np.array(radius_holder)[:,1]))
# plt.title("Correlation Lenght Fit With Exponential Only")
# plt.xlabel("Sample Frame")
# plt.ylabel("Correlation Lenght (units of lattice)")
# plt.show()


# plt.plot((np.array(popt1)[:,0]))
# plt.title("Correlation Lenght Fit With Exponential Only")
# plt.xlabel("Sample Frame")
# plt.ylabel("Correlation Lenght (units of lattice)")
# plt.show()





##############################################################################################
####correlation length stuff



# pc_holder=[]
# radius_holder=[]
# radius_holder2=[]
# radi=[]
# popt1=[]

# for i in range(len(VA_data_type_O)):
#    temp1, temp2, temp3 = paircorrelation3D_a(x_per_frame[i],y_per_frame[i],z_per_frame[i],sizeN/2.,CM_ar[i],c_per_frame[i],dr=0.5)
#    pc_holder.append(temp1)
#    radius_holder.append(temp2)
#    radi.append(temp3)

'''
#radius_holder2.append(temp3)
for i in range(len(pc_holder)):
   plt.plot(radi[i],pc_holder[i],'ro')
   plt.plot(radi[i],epo1(radi[i],radius_holder[i][0],radius_holder[i][1]))
   plt.title("Correlation Function")
   plt.xlabel("units")
   plt.ylabel("Averaged Correlation Function")
   plt.yscale("log")
   #plt.xscale("log")
   plt.show()

'''
# expo_correlation = 1./(np.array(radius_holder)[:,1])
# plt.plot(1./(np.array(radius_holder)[:,1]))
# plt.title("Correlation Length w Exponential Over Samples")
# plt.xlabel("Sample #")
# plt.ylabel("Correlation Length (lattice units)")
# plt.show()






# for i in range(len(distance)):
#    plt.hist(pc_holder[i])

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

# animate_DT = []

# for j in range(0,samples[0]):
#     animate_DT.append(np.reshape(VA_data_type_O[j],(T_P,4)))













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
    popt, pcov = curve_fit(fit_MSD,list(range(1,len(d)+1))[:10],d[:10],p0=[1,1],maxfev=10000)
    #popt, pcov = curve_fit(fit_MSD,range(1,len(d)+1)[:int(len(range(1,len(d)+1))*thresh)],d[:int(len(range(1,len(d)+1))*thresh)],p0=[1,1],maxfev=10000)
    return [popt,pcov]

def MSD_Calculation(x):
    track_ptype.append(x)

    #msd
    msd_t_mean = np.zeros(len(x))


    fit_verbose = []
    msd_s_mean = []
    msd_s_mean_v = []


    for i in range(len(x)):
        msd_t_mean_p = np.zeros(len(x[i]))
        msd_dec = []
        
        for j in range(len(x[i])):
            pp = x[i][j]
            
            msd_t_mean_p[j] = MSD_tavg(pp[:,0],pp[:,1],pp[:,2],1,sizeN)
            tttt=track_decomp(pp[:,0],pp[:,1],pp[:,2],1,sizeN)
            msd_dec.append(tttt)
        


        msd_dec = np.mean(np.asarray(msd_dec),axis = 0)
        msd_s_mean.append(msd_dec)
                        
        msd_t_mean[i] = np.mean(msd_t_mean_p)

    msd_s_meana_cm.append(np.mean(np.asarray(msd_s_mean),axis = 0))

    fit_msdaa_cm.append(fit_decom(np.mean(np.asarray(msd_s_mean),axis = 0))[0])

    msd_t_meana_cm[k] = np.mean(msd_t_mean)
    msd_t_mean_verba_cm.append(msd_t_mean)
    return 0


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


msd_t_meana_cm = np.zeros(len(sizeM))
msd_t_mean_verba_cm = []
    

fit_msdaa_cm = []


msd_s_meana_cm = []






#global holder array for angle stuff
ang_global = []


for k in range(len(sizeM)):


    #structured arrays for msd calculations
    n_ordered_arr = np.zeros((sizeM[k],sizeP[k],samples[0],4))
    cm_ordered_arr = np.zeros((1,sizeP[k],samples[0],3))


    #structured arrays for angle between line n trajectory
    #setup as n_ordered_arr but with size ((sizeM[k],sizeP[k],samples[0]-2))

    ang_ordered_arr = np.zeros((sizeM[k],sizeP[k],samples[0]-2))


    for i in range(len(VA_data_type_O)):
        for j in range(len(VA_data_type_O[i][k])):
            for l in range(len(VA_data_type_O[i][k][j])):
                n_ordered_arr[l][j][i] = VA_data_type_O[i][k][j][l]
            cm_ordered_arr[0][j][i] = cm(VA_data_type_O[i][k][j][:,0],VA_data_type_O[i][k][j][:,1],VA_data_type_O[i][k][j][:,2],sizeN)




    #do the angle calculations
    for i in range(len(ang_ordered_arr)):
        for j in range(len(ang_ordered_arr[i])):
            temp_array = np.asarray(n_ordered_arr[i][j])
            ang_ordered_arr[i][j] = angle_trajectory_3d(temp_array[:,0],temp_array[:,1],temp_array[:,2])

    ang_global.append(ang_ordered_arr)  
    # outdated code to save data to file for angle -> now do it in this script: TODO -> make a modular way for this to occur
    # import csv
    # with open('employee_file_{0}.csv'.format(k), mode='w') as employee_file:
    #     employee_writer = csv.writer(employee_file, delimiter=',')

    #     for j in n_ordered_arr[0]:
    #         for i in n_ordered_arr[0][0]:
    #             employee_writer.writerow(i)




    MSD_Calculation(cm_ordered_arr)
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

    fit_msdaa.append(fit_decom(np.mean(np.asarray(msd_s_mean),axis = 0))[k])

    msd_t_meana[k] = np.mean(msd_t_mean)
    msd_t_mean_verba.append(msd_t_mean)




#plotting angle stuff


for i in range(len(ang_global[0])):
    for j in range(len(ang_global[0][i])):
        hist_polar(ang_global[0][i][j])



hist_polar(ang_global[1][0][0])











def msd_line(D_A,length):
    return D_A[0]*list(range(length))**D_A[1]

def msd_line2(D_A,length):
    return D_A[0]*np.array((range(length)))**D_A[1]

def msd_plot(m,c12,a,b,verbose = False,cm = False): #a,b = expanded values of fit_msdaa
    #Verbose: all msd or averaged ; True, False
    if verbose:
        for i in range(len(m)):
            for j in m[i]:
                plt.plot(list(range(1,len(j)+1)),j)
            plt.ylabel("MSD")
            plt.xlabel("Tau")
            plt.title("P = {0}, M = {1}".format(sizeP[i],sizeM[i]))
            plt.yscale("log")
            plt.xscale("log")
            plt.show()
    else:
        for i in range(len(m)):
            plt.plot(list(range(1,len(m[i])+1)),m[i],label="{0}".format(i+1))
            plt.plot(list(range(1,len(c12[i])+1)),c12[i])
            if cm:
                plt.plot(msd_line(fit_msdaa_cm[0],len(m[i])))  
            else:
                plt.plot(msd_line(fit_msdaa[0],len(m[i])))
                #plt.plot(msd_line2([a,b],len(m[i])))
            #plt.plot(msd_line(fit_msdaa_cm[0],len(m[i])))
        plt.ylabel("Log10 MSD (units of simulation box)")
        plt.xlabel("Log10 Tau (timestep of trajectory)")
        #plt.title("P = {0}, M = {1}".format(sizeP,sizeM))
        plt.yscale("log")
        plt.xscale("log")
        #plt.legend()
        plt.savefig("MSD curve")
        plt.show()
    return


# fig = plt.figure()
# ax = fig.add_subplot(211)
# ax2 = fig.add_subplot(212)

# for i in range(len(msd_s_meana)):
#     ax.plot(list(range(1,len(msd_s_meana[i])+1)),msd_s_meana[i],label="{0}".format(i+1))
#     ax.plot(list(range(1,len(msd_s_meana_cm[i])+1)),msd_s_meana_cm[i])

#     ax.plot(msd_line(fit_msdaa[0],len(msd_s_meana[i])))
#         #plt.plot(msd_line2([a,b],len(m[i])))
#     #plt.plot(msd_line(fit_msdaa_cm[0],len(m[i])))
# ax2.plot(fit_verbosed[0])
# ax.annotate('A', xy=(-40, -10 + ax2.bbox.height), xycoords="axes pixels", fontsize=25, weight = 'bold')
# ax2.annotate('B', xy=(-40, -10 + ax2.bbox.height), xycoords="axes pixels", fontsize=25, weight = 'bold')
# ax.set_ylabel("Log10 MSD (au)")
# ax.set_xlabel("Log10 Tau")
# #plt.title("P = {0}, M = {1}".format(sizeP,sizeM))
# ax.set_yscale("log")
# ax.set_xscale("log")

# ax2.set_ylabel("Diffusion Coefficient (au)")
# ax2.set_xlabel("Monomer Index")
# #plt.title("P = {0}, M = {1}".format(sizeP,sizeM))
# # ax2.yscale("log")
# # ax2.xscale("log")
# #plt.legend()
# fig.tight_layout()
# fig.savefig("time")
# fig.show()



def create_box_plot(box_data,tick_list,y_label = "",x_label = "",y_lim = (),title = ""):
    ticks = tick_list
    plt.boxplot(box_data,positions = list(range(1,len(tick_list)+1)))
    for i in range(1,len(tick_list)+1):
        y = box_data[i-1]
        x = np.random.normal(i, 0.04, size=len(y))
        plt.plot(x, y, 'r.', alpha=0.2)
    try:
        plt.ylim(y_lim)
    except:
        print("Warning: y_lim not valid")
    plt.xticks(range(1, len(ticks) * 1 + 1, 1), ticks)
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
        print("Fitted Mean: {0} +/- {1}".format(gmix.means_[:,1],np.sqrt(gmix.covars_[:,1])))
        print("Fitted Mean(normal): {0} +/- {1}".format(np.exp(gmix.means_[:,1]),np.exp(gmix.means_[:,1])*np.sqrt(gmix.covars_[:,1])))
    else:
        print("Fitted Mean: {0} +/- {1}".format(gmix.means_[:,1],np.sqrt(gmix.covars_[:,1])))
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
        print("Warning: x_limit is invalid")
    plt.show()

    return




def plot_distance_time(distances):

    fig, ax = plt.subplots()

    ax.plot(np.arange(0,R,(R/(samples[0]-1))),distances[1:])
    ax.set_title("Total distance between all monomers over time")
    ax.set_xlabel("Simulation time (arbritrary units)")
    ax.set_ylabel("Total distance (units of Simulation Box)")
    fig.savefig("Distance over time")
    plt.tight_layout()
    plt.show()

    return [np.arange(0,R,(R/(samples[0]-1))),distances[1:]]

def plot_correlation_time(correlation):

    fig, ax = plt.subplots()

    ax.plot(np.arange(0,R,(R/(samples[0]-1))),correlation[1:])
    ax.set_title("Correlation length of system over time")
    ax.set_xlabel("Simulation time (arbritrary units)")
    ax.set_ylabel("Correlation length (units of simulation box)")
    fig.savefig("Correlation over time")
    plt.tight_layout()
    plt.show()

    return [np.arange(0,R,(R/(samples[0]-1))),correlation[1:]]




# return_dis = plot_distance_time(distance_arr)
# return_corr = plot_correlation_time(expo_correlation)








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






















def normalize_colors(color_float, max_c = 10., min_c = -10., range_min = 0., range_max = 1.):

    return (range_max - range_min) * ((color_float - min_c))/(max_c - min_c) + range_min


import matplotlib.cm as cm_color
def hi(val = 10):
    #############################################################################################
    #
    ##plotting and animations
    #
    #############################################################################################

    #Set up formatting for the movie files
    Writer = matplotlib.animation.writers['ffmpeg']
    writer = Writer(fps=10, metadata=dict(artist='Baljyot Parmar, all rights reserved'), bitrate=1800)

    norm = colors.Normalize(vmin=-10, vmax=10)
    m = cm_color.ScalarMappable(norm=norm, cmap="winter")


    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    title = ax.set_title('3D Test')


    def update_graph(num):
        i = animate_DT[num]

        graph._offsets3d = (i[:,0],i[:,1],i[:,2])


        title.set_text('3D Test, Sample={}'.format(num))

    graph = ax.scatter(animate_DT[0][:,0],animate_DT[0][:,1],animate_DT[0][:,2]) # for color coded with charge add c=m.to_rgba(animate_DT[0][:,3]) to plot argument

    ax.set_xbound(lower=-1,upper=sizeN+1)
    ax.set_ybound(lower=-1,upper=sizeN+1)
    ax.set_zbound(lower=-1,upper=sizeN+1)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")
    ani = matplotlib.animation.FuncAnimation(fig, update_graph, list(range(1,samples[0])), blit=False)
    ##ani.save('/Users/baljyot/Documents/Polymer_Output/new_ani.mp4',writer=writer)
    ani.save('{0}/new_ani.mp4'.format(os.getcwd()),writer=writer)
    plt.show()




    fig = plt.figure()
    ax = fig.add_subplot(221,projection = '3d')
    ax.set_xlim([0, sizeN])
    ax.set_ylim([0, sizeN])
    ax.set_zlim([0, sizeN])

    ax2 = fig.add_subplot(222,projection = '3d')
    ax2.set_xlim([0, sizeN])
    ax2.set_ylim([0, sizeN])
    ax2.set_zlim([0, sizeN])

    ax3 = fig.add_subplot(223)
    ax4 = fig.add_subplot(224)

    ax.scatter(animate_DT[0][:,0],animate_DT[0][:,1],animate_DT[0][:,2]) #color add , c=m.to_rgba(animate_DT[0][:,3])
    ax.set_title("Simulation time = 0")
    ax2.scatter(animate_DT[samples[0]-1][:,0],animate_DT[samples[0]-1][:,1],animate_DT[samples[0]-1][:,2]) # for color coded with charge add c=m.to_rgba(animate_DT[samples[0]-1][:,3]) to plot argument
    ax2.set_title("Simulation time = {0}".format(R))
    ax3.plot(np.arange(0,R,(R/(samples[0]-1))),distance_arr[1:])
    ax3.set_xticks(np.linspace(0,R,3))
    ax4.plot(np.arange(0,R,(R/(samples[0]-1))),1./(np.array(radius_holder)[:,1])[1:])
    ax4.set_xticks(np.linspace(0,R,3))
    ax3.set_xlabel("Simulation time (au)")
    ax3.set_ylabel("Total distance (au)")
    ax4.set_xlabel("Simulation time (au)")
    ax4.set_ylabel("Correlation length (au)")
    ax.annotate('A', xy=(-40, val + ax2.bbox.height), xycoords="axes pixels", fontsize=25, weight = 'bold')
    ax2.annotate('B', xy=(-40, val + ax2.bbox.height), xycoords="axes pixels", fontsize=25, weight = 'bold')
    ax3.annotate('C', xy=(-40, val + ax2.bbox.height), xycoords="axes pixels", fontsize=25, weight = 'bold')
    ax4.annotate('D', xy=(-40, val + ax2.bbox.height), xycoords="axes pixels", fontsize=25, weight = 'bold')
    plt.tight_layout()
    plt.savefig("Panel")

    plt.plot(animate_DT[0][:,0],animate_DT[0][:,1],'ro')
    plt.show()
    '''
    '''
    for k in range(0,samples[0]):
        for i in VA_data_type_O[k]:
            for j in i:
                ax.scatter(j[:,0],j[:,1],j[:,2]) #for color add ,c=j[:,3] to argument
        #ax.plot(j[:,0],j[:,1],j[:,2],'-b')
        ax.set_xbound(lower=0,upper=sizeN)
        ax.set_ybound(lower=0,upper=sizeN)
        ax.set_zbound(lower=0,upper=sizeN)
    plt.show()

    return 





