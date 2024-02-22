"""
Created on Mon Dec  5 11:58:04 2022

@author: c.d.hamilton@lancaster.ac.uk

read results and perform post-processing analysis
"""

import numpy as np

from matplotlib import pyplot as plt

# to use Latex font type
plt.rcParams['text.usetex'] = True

# load npy file
data = np.load("SMVEM_1000s_method3.npy", allow_pickle=True)

niter = len(data)

time = []
sun = []
mercury = []
venus = []
earth = []
mars = []

# create lists with time and energy information only from data
for it in range(0,niter):
    time.append(data[it][0])
    sun.append(data[it][1])
    mercury.append(data[it][2])
    venus.append(data[it][3])
    earth.append(data[it][4])
    mars.append(data[it][5])
    
    sunx = []
    suny = []
    mercuryx = []
    mercuryy = []
    venusx = []
    venusy = []
    earthx = []
    earthy = []
    marsx = []
    marsy = []
    PE = []
    KE = []
    Tot_E = []
    LM = []
    AM = []
    Tot_EF = []
    LMF = []
    AMF = []
    
    for i in range(0, len(mercury)):
        sunx.append(sun[i].position[0])
        suny.append(sun[i].position[1])
        mercuryx.append(mercury[i].position[0])
        mercuryy.append(mercury[i].position[1])
        venusx.append(venus[i].position[0])
        venusy.append(venus[i].position[1])
        earthx.append(earth[i].position[0])
        earthy.append(earth[i].position[1])
        marsx.append(mars[i].position[0])
        marsy.append(mars[i].position[1])
        PE.append(sun[i].PE + mercury[i].PE + venus[i].PE + earth[i].PE + mars[i].PE)
        KE.append(sun[i].KE + mercury[i].KE + venus[i].KE + earth[i].KE + mars[i].KE)
        Tot_E.append((((0.5 * sun[i].PE ) + sun[i].KE) + ((0.5 * mercury[i].PE ) + mercury[i].KE) + ((0.5 * venus[i].PE ) + venus[i].KE) + ((0.5 * earth[i].PE ) + earth[i].KE) + ((0.5 * mars[i].PE ) + mars[i].KE)))
        LM.append(np.linalg.norm(sun[i].LM) + np.linalg.norm(mercury[i].LM) + np.linalg.norm(venus[i].LM) + np.linalg.norm(earth[i].LM) + np.linalg.norm(mars[i].LM))
        AM.append(np.linalg.norm(sun[i].AM) + np.linalg.norm(mercury[i].AM) + np.linalg.norm(venus[i].AM) + np.linalg.norm(earth[i].AM) + np.linalg.norm(mars[i].AM))
        Tot_EF.append((((0.5 * sun[i].PE ) + sun[i].KE) + ((0.5 * mercury[i].PE ) + mercury[i].KE) + ((0.5 * venus[i].PE ) + venus[i].KE) + ((0.5 * earth[i].PE ) + earth[i].KE) + ((0.5 * mars[i].PE ) + mars[i].KE)) / (((0.5 * sun[0].PE ) + sun[0].KE) + ((0.5 * mercury[0].PE ) + mercury[0].KE) + ((0.5 * venus[0].PE ) + venus[0].KE) + ((0.5 * earth[0].PE ) + earth[0].KE) + ((0.5 * mars[0].PE ) + mars[0].KE)))
        LMF.append((np.linalg.norm(sun[i].LM) + np.linalg.norm(mercury[i].LM) + np.linalg.norm(venus[i].LM) + np.linalg.norm(earth[i].LM) + np.linalg.norm(mars[i].LM)) / (np.linalg.norm(sun[0].LM) + np.linalg.norm(mercury[0].LM) + np.linalg.norm(venus[0].LM) + np.linalg.norm(earth[0].LM) + np.linalg.norm(mars[0].LM)))
        AMF.append((np.linalg.norm(sun[i].AM) + np.linalg.norm(mercury[i].AM) + np.linalg.norm(venus[i].AM) + np.linalg.norm(earth[i].AM) + np.linalg.norm(mars[i].AM)) / (np.linalg.norm(sun[0].AM) + np.linalg.norm(mercury[0].AM) + np.linalg.norm(venus[0].AM) + np.linalg.norm(earth[0].AM) + np.linalg.norm(mars[0].AM)))
        
#Plot of the planetary positions over the course of the simulation
"""
fig = plt.figure(figsize=(10,10), dpi = 250)
plt.rc('font', size = 20)
plt.plot(sunx, suny, 'y', marker = "x", label = 'Sun')
plt.plot(mercuryx, mercuryy, 'k', label = 'Mercury')
plt.plot(venusx, venusy, 'g', label = 'Venus')
plt.plot(earthx, earthy, 'b', label = 'Earth')
plt.plot(marsx, marsy, 'r', label = 'Mars')
plt.xlabel("x [m]", fontsize=20)
plt.ylabel("y [m]", fontsize=20)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.legend(loc = 'lower right', fontsize='large')
plt.savefig('Orbits_2000s_method3.jpeg')
plt.show()
"""

#Plot of ratios of total energy, linear momentum and angular momentum
fig = plt.figure(figsize=(16,8))
# first sublot will be time vs total energy
ax1 = fig.add_subplot(3,1,1)
ax1.set_ylabel("TE",fontsize=20)
ax1.tick_params(labelsize=20)
ax1.ticklabel_format(useOffset=False)
plot1, = ax1.plot(time , Tot_EF,'r')
# second sublot will be time vs linear momentum
ax2 = fig.add_subplot(3,1,2)
ax2.set_ylabel("LM",fontsize=20)
ax2.tick_params(labelsize=20)
ax2.ticklabel_format(useOffset=False)
plot2, = ax2.plot(time, LMF,'g')
# third sublot will be time vs angular momentum
ax3 = fig.add_subplot(3,1,3)
ax3.set_xlabel("Time [s]",fontsize=20)
ax3.set_ylabel("AM",fontsize=20)
ax3.tick_params(labelsize=20)
ax3.ticklabel_format(useOffset=False)
plot3, = ax3.plot(time, AMF,'b.')
plt.savefig('Properties_2000s_method3.jpeg')
