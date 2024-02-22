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
data = np.load("SM_1000s_7600_method2.npy", allow_pickle=True)

niter = len(data)

time = []
sun = []
mercury = []

# create lists with time and energy information only from data
for it in range(0,niter):
    time.append(data[it][0])
    sun.append(data[it][1])
    mercury.append(data[it][2])
    
    sunx = []
    suny = []
    mercuryx = []
    mercuryy = []
    
    
    for i in range(0, len(mercury)):
        sunx.append(sun[i].position[0])
        suny.append(sun[i].position[1])
        mercuryx.append(mercury[i].position[0])
        mercuryy.append(mercury[i].position[1])
        


fig = plt.figure(figsize=(10,10), dpi = 250)
plt.rc('font', size = 20)
plt.plot(sunx, suny, 'y', marker = "x", label = 'Sun')
plt.plot(mercuryx, mercuryy, 'k', label = 'Mercury')
plt.plot(mercury[0].position[0], mercury[0].position[1], marker = '.', label = 'Mercury Start', markersize = 50, mec = 'b')
plt.plot(mercury[1519].position[0], mercury[1519].position[1], marker = '.', label = 'Mercury End', markersize = 50, mec = 'r')
plt.xlabel("x [m]", fontsize=20) 
plt.ylabel("y [m]", fontsize=20)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.legend(loc = 'upper right', fontsize='large')
plt.savefig('TwoBody_1000s_7600_method2.jpeg')
plt.show()