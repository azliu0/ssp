import numpy as np
from math import sin, cos

import ephem

obs = ephem.Observer() 
obs.lon = "-81:24:52.9" #longitude of obs. 
obs.lat = "36:15:09.7" #lat of obs. 
obs.elev = 922. #elevation of obs, in meters (not a string!) 
obs.date = "2020/07/05 07:08:00" #(UTC date/time of observation) 

line = "2003GE42,e,27.8543085,165.6469288,101.353156,2.63408,,0.375418379,8.4377868,07/05.29722/2020,2000,,"

'''
Field1 = name
Field2 = e (for Elliptical body)
Field3 = inclination, degrees
Field4 = ascending node, degrees
Field5 = argument of perihelion, degrees
Field6 = semi-major axis, au
Field7 = blank (mean motion, n, optional; note the extra comma)
Field8 = eccentricity
Field9 = mean anomaly, degrees
Field10 = epoch for mean anomaly (MM/DD.ddd/YYYY - turn UTC time into adecimal day)
Field11 = 2000.0 (the equinox)
'''

asteroid = ephem.readdb(line)
asteroid.compute(obs)
print(asteroid.a_ra,asteroid.a_dec)