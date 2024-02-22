"""
Created on Mon Dec  5 11:57:00 2022

@author: c.d.hamilton@lancaster.ac.uk

creates instances of the Particle class, advances the system in time and saves results for use in two-body solution
"""

import numpy as np
from particle import Particles
import copy
from astropy.time import Time
from astropy.coordinates import get_body_barycentric_posvel
from poliastro import constants
from astropy.constants import G
from spiceypy import sxform, mxvg

# get the time at 5pm on 5th December 2022
t = Time("2022-12-05 17:00:00.0", scale="tdb")

# get transformation matrix to the ecliptic (use time in Julian Days)
trans = sxform("J2000", "ECLIPJ2000", t.jd)

## Sun values
sunpos, sunvel = get_body_barycentric_posvel("sun", t, ephemeris="jpl")
# make a "state vector" of positions and velocities (in metres and metres/second, respectively)
sunstatevec = [
    sunpos.xyz[0].to("m").value,
    sunpos.xyz[1].to("m").value,
    sunpos.xyz[2].to("m").value,
    sunvel.xyz[0].to("m/s").value,
    sunvel.xyz[1].to("m/s").value,
    sunvel.xyz[2].to("m/s").value,]
# transform state vector to ecliptic
sunstatevececl = mxvg(trans, sunstatevec)
# get positions and velocities
sunposition = [sunstatevececl[0], sunstatevececl[1], sunstatevececl[2]]
sunvelocity = [sunstatevececl[3], sunstatevececl[4], sunstatevececl[5]]

## Mercury values
mercurypos, mercuryvel = get_body_barycentric_posvel("mercury", t, ephemeris="jpl")
# make a "state vector" of positions and velocities (in metres and metres/second, respectively)
mercurystatevec = [
    mercurypos.xyz[0].to("m").value,
    mercurypos.xyz[1].to("m").value,
    mercurypos.xyz[2].to("m").value,
    mercuryvel.xyz[0].to("m/s").value,
    mercuryvel.xyz[1].to("m/s").value,
    mercuryvel.xyz[2].to("m/s").value,]
# transform state vector to ecliptic
mercurystatevececl = mxvg(trans, mercurystatevec)
# get positions and velocities
mercuryposition = [mercurystatevececl[0], mercurystatevececl[1], mercurystatevececl[2]]
mercuryvelocity = [mercurystatevececl[3], mercurystatevececl[4], mercurystatevececl[5]]

Sun = Particles(position = np.array(sunposition, dtype=float),
                velocity = np.array(sunvelocity, dtype=float),
                acceleration = np.array([0, 0, 0], dtype=float),
                name = "Sun",
                mass = (constants.GM_sun / G).value,
                KE = 0.5 * (constants.GM_sun / G).value * (np.linalg.norm(sunvelocity) ** 2),
                PE = 0,
                LM = 0,
                AM = 0,
                LA = np.array([0, 0, 0], dtype=float))


Mercury = Particles(position = np.array(mercuryposition, dtype=float),
                velocity = np.array(mercuryvelocity, dtype=float),
                acceleration = np.array([0, 0, 0], dtype=float),
                name = "Mercury",
                mass = (constants.GM_mercury / G).value,
                KE = 0.5 * (constants.GM_mercury / G).value * (np.linalg.norm(mercuryvelocity) ** 2),
                PE = 0,
                LM = 0,
                AM = 0,
                LA = np.array([0, 0, 0], dtype=float))

bodies = [Sun, Mercury] #list of all bodies to be included in the simulation
i=1
Data = []

for i in range(0, 760):
    for body in bodies:
        body.updateGravitationalAcceleration(bodies)
        body.update(10000, 3)
    time = i * 10000
    if (i%760) == 0:
        print(i/760)
    if i == 0 or i % 5 == 0:
        Data.append([time, copy.deepcopy(Sun), copy.deepcopy(Mercury)])
        i += 1
    
np.save("SM_10000s_760", Data, allow_pickle=True)