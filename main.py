"""
Created on Mon Dec  5 11:57:00 2022

@author: c.d.hamilton@lancaster.ac.uk

creates instances of the Particle class, advances the system in time and saves results
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

## Venus values
venuspos, venusvel = get_body_barycentric_posvel("venus", t, ephemeris="jpl")
# make a "state vector" of positions and velocities (in metres and metres/second, respectively)
venusstatevec = [
    venuspos.xyz[0].to("m").value,
    venuspos.xyz[1].to("m").value,
    venuspos.xyz[2].to("m").value,
    venusvel.xyz[0].to("m/s").value,
    venusvel.xyz[1].to("m/s").value,
    venusvel.xyz[2].to("m/s").value,]
# transform state vector to ecliptic
venusstatevececl = mxvg(trans, venusstatevec)
# get positions and velocities
venusposition = [venusstatevececl[0], venusstatevececl[1], venusstatevececl[2]]
venusvelocity = [venusstatevececl[3], venusstatevececl[4], venusstatevececl[5]]

## Earth values
earthpos, earthvel = get_body_barycentric_posvel("earth", t, ephemeris="jpl")
# make a "state vector" of positions and velocities (in metres and metres/second, respectively)
earthstatevec = [
    earthpos.xyz[0].to("m").value,
    earthpos.xyz[1].to("m").value,
    earthpos.xyz[2].to("m").value,
    earthvel.xyz[0].to("m/s").value,
    earthvel.xyz[1].to("m/s").value,
    earthvel.xyz[2].to("m/s").value,]
# transform state vector to ecliptic
earthstatevececl = mxvg(trans, earthstatevec)
# get positions and velocities
earthposition = [earthstatevececl[0], earthstatevececl[1], earthstatevececl[2]]
earthvelocity = [earthstatevececl[3], earthstatevececl[4], earthstatevececl[5]]

## Mars values
marspos, marsvel = get_body_barycentric_posvel("mars", t, ephemeris="jpl")
# make a "state vector" of positions and velocities (in metres and metres/second, respectively)
marsstatevec = [
    marspos.xyz[0].to("m").value,
    marspos.xyz[1].to("m").value,
    marspos.xyz[2].to("m").value,
    marsvel.xyz[0].to("m/s").value,
    marsvel.xyz[1].to("m/s").value,
    marsvel.xyz[2].to("m/s").value,]
# transform state vector to ecliptic
marsstatevececl = mxvg(trans, marsstatevec)
# get positions and velocities
marsposition = [marsstatevececl[0], marsstatevececl[1], marsstatevececl[2]]
marsvelocity = [marsstatevececl[3], marsstatevececl[4], marsstatevececl[5]]

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

Venus = Particles(position = np.array(venusposition, dtype=float),
                velocity = np.array(venusvelocity, dtype=float),
                acceleration = np.array([0, 0, 0], dtype=float),
                name = "Venus",
                mass = (constants.GM_venus / G).value,
                KE = 0.5 * (constants.GM_venus / G).value * (np.linalg.norm(venusvelocity) ** 2),
                PE = 0,
                LM = 0,
                AM = 0,
                LA = np.array([0, 0, 0], dtype=float))

Earth = Particles(position = np.array(earthposition, dtype=float),
                velocity = np.array(earthvelocity, dtype=float),
                acceleration = np.array([0, 0, 0], dtype=float),
                name = "Earth",
                mass = (constants.GM_earth / G).value,
                KE = 0.5 * (constants.GM_earth / G).value * (np.linalg.norm(earthvelocity) ** 2),
                PE = 0,
                LM = 0,
                AM = 0,
                LA = np.array([0, 0, 0], dtype=float))

Mars = Particles(position = np.array(marsposition, dtype=float),
                velocity = np.array(marsvelocity, dtype=float),
                acceleration = np.array([0, 0, 0], dtype=float),
                name = "Mars",
                mass = (constants.GM_mars / G).value,
                KE = 0.5 * (constants.GM_mars / G).value * (np.linalg.norm(marsvelocity) ** 2),
                PE = 0,
                LM = 0,
                AM = 0,
                LA = np.array([0, 0, 0], dtype=float))

bodies = [Sun, Mercury, Venus, Earth, Mars] #list of all bodies to be included in the simulation
i=1
Data = []

for i in range(0, 30000):
    for body in bodies:
        body.updateGravitationalAcceleration(bodies)
        body.update(2000, 3)
        body.updateKE()
        body.updatePE(bodies)
        body.updateLM()
        body.updateAM()
    time = i * 2000
    if (i%300) == 0:
        print(i/300)
    if i == 0 or i % 1000 == 0:
        Data.append([time, copy.deepcopy(Sun), copy.deepcopy(Mercury), copy.deepcopy(Venus), copy.deepcopy(Earth), copy.deepcopy(Mars)])
        i += 1
    
np.save("SMVEM_2000s_method3", Data, allow_pickle=True)