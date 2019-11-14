### adapted from https://github.com/cfinch/Shocksolution_Examples/tree/master/PairCorrelation

import numpy as np
def epo(x,*p0):
    return (1./(x**p0[0]))*p0[1]*np.exp(-x*p0[2]) #+ p0[3]
def epo1(x,*p0):
    return p0[0]*np.exp(-x*p0[1])

def dist(x,y,z,c,N):
    t1=np.minimum(abs(c[0]-x),abs(c[0]-x - N))
    tx=np.minimum(t1,abs(c[0]-x + N))
    t2=np.minimum(abs(c[1]-y),abs(c[1]-y - N))
    ty=np.minimum(t2,abs(c[1]-y + N))
    t3=np.minimum(abs(c[2]-z),abs(c[2]-z - N))
    tz=np.minimum(t3,abs(c[2]-z + N))
    temp=np.sqrt((tx)**2 + (ty)**2 + (tz)**2)
    return temp  
      
def paircorrelation3D(x,y,z,S,CM,chg,dr=0.1):
    from numpy import zeros, sqrt, pi, mean, arange, histogram, min, digitize, sum, fft
    from scipy.optimize import curve_fit
    import pylab as py
    
    #rmax=min([min([CM[0],(S-CM)[0]]),min([CM[1],(S-CM)[1]])])  #this is not taking into account periodic boundary conditions and only looks at a set window
    #to account for periodic boundary we allow rmax to be of the lattice size from the reference of the center of mass
    
    rmax=S

    d= dist(x,y,z,CM,S)
    edges=arange(0.,rmax + 1.1*dr, dr)
    (results,bins)=histogram(d,bins=edges,normed=False)
    temp=digitize(d,bins,right=True)
    num_increments = len(edges) - 1
    radii = zeros(num_increments)
    numberDensity = float(len(x)) / S**3
    g = results/numberDensity
    g_average = zeros(num_increments)
    avg_chg = zeros(num_increments)
    gg=zeros(num_increments)
    for i in range(num_increments):
        avg_chg[i]=sum(chg[i==temp])
        radii[i] = (edges[i] + edges[i+1]) / 2.
        rOuter = edges[i + 1]
        rInner = edges[i]
        g_average[i] = mean(g[i]) / ((4.*pi * (rOuter**3 - rInner**3))/3.)
        gg[i]=g[i]/((4.*pi * (rOuter**3 - rInner**3))/3.)
    max_v=pi*rmax**2
    k=0!=g_average
    kkk=g_average[k]
    kk=int(round(0.1*len(kkk)))
    iq=max_v*fft.fft(g_average[k])
    iqf=fft.fftfreq(len(g_average[k]),dr)

    popt1 = [1,1]
    pcov1 = [[1,1],[1,1]]
    #popt, pcov = curve_fit(epo,radii[k],g_average[k],p0=[1,1,1],maxfev=100000)
    try:
        popt1, pcov1 = curve_fit(epo1,radii[k][int(round(len(radii[k])*0.1)):],g_average[k][int(round(len(radii[k])*0.1)):],p0=[1,1],maxfev=100000)
    except:
        pass
    #print popt[0],popt[1],1./popt[2], sqrt(pcov[0,0]),sqrt(pcov[1,1]),1./sqrt(pcov[2,2])
    #print popt1[0], popt1[1], 1./popt1[1], sqrt(pcov1[0,0]),sqrt(pcov1[1,1])
    #return (g_average, radii,rmax,avg_chg,iq,iqf,popt1,pcov1)
    
    return g_average,popt1,gg,radii

def paircorrelation3D_a(x,y,z,S,CM,chg,dr=0.1):
    from numpy import zeros, sqrt, pi, mean, arange, histogram, min, digitize, sum, fft
    from scipy.optimize import curve_fit
    import pylab as py
    
    #rmax=min([min([CM[0],(S-CM)[0]]),min([CM[1],(S-CM)[1]])])  #this is not taking into account periodic boundary conditions and only looks at a set window
    #to account for periodic boundary we allow rmax to be of the lattice size from the reference of the center of mass
    
    rmax=S
    edges=arange(0.,rmax + 1.1*dr, dr)
    num_increments = len(edges) - 1

    g = zeros([len(x), num_increments])
    g1 = zeros([len(x), num_increments])
    numberDensity = float(len(x)) / S**3
    for p in range(len(x)):
        d= dist(x,y,z,[x[p],y[p],z[p]],S)
        (results,bins)=histogram(d,bins=edges,normed=False)
        g[p,:] = results/numberDensity
        g1[p,:] = results/numberDensity
    
    radii = zeros(num_increments)
    g_average = zeros(num_increments)

    for i in range(num_increments):
        radii[i] = (edges[i] + edges[i+1]) / 2.
        rOuter = edges[i + 1]
        rInner = edges[i]
        g_average[i] = mean(g[:,i]) / ((4.*pi * (rOuter**3 - rInner**3))/3.)
        g1[:,i] = g1[:,i]/((4.*pi * (rOuter**3 - rInner**3))/3.)
    max_v=pi*rmax**2
    k=0!=g_average
    kkk=g_average[k]
    kk=int(round(0.1*len(kkk)))
    iq=max_v*fft.fft(g_average[k])
    iqf=fft.fftfreq(len(g_average[k]),dr)
#popt, pcov = curve_fit(epo,radii[k],g_average[k],p0=[1,1,1],maxfev=100000)
    try:
        popt1, pcov1 = curve_fit(epo1,radii[k][:int(round(len(radii[k])*0.5))],g_average[k][:int(round(len(radii[k])*0.5))],p0=[1,1],maxfev=100000)
    except:
        popt1 = [0,0,0]
    #print popt[0],popt[1],1./popt[2], sqrt(pcov[0,0]),sqrt(pcov[1,1]),1./sqrt(pcov[2,2])
    #print popt1[0], popt1[1], 1./popt1[1], sqrt(pcov1[0,0]),sqrt(pcov1[1,1])
    #return (g_average, radii,rmax,avg_chg,iq,iqf,popt1,pcov1)
    
    return g_average,popt1,radii#,popt
