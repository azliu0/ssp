from vpython import vector, sphere, curve, rate, color
from math import radians, sin, cos, sqrt, pi
import numpy as np

#initializing orbital elements
a = 2.773017979589484
e = 0.1750074901308245
M = radians(336.0050001501443)
Oprime = radians(108.032597191534)
iprime = radians(16.34548466739393)
wprime = radians(74.95130563682554)


#newton-gauss method to solve kepler equation
def solvekep(M):
    Eguess = M
    Mguess = Eguess - e*sin(Eguess)
    Mtrue = M
    while (abs(Mguess - Mtrue) > 1e-004):
        Mguess = Eguess - e*sin(Eguess)
        Eguess = Eguess - (Eguess - e*sin(Eguess) - Mtrue) / (1 - e*cos(Eguess))
    return Eguess


#period
sqrtmu = 0.01720209895
mu = sqrtmu**2
time = 0
period = sqrt(4*pi**2*a**3/mu)

#ecliptic coordinates
r1ecliptic = np.array([0, 0, 0])
Mtrue = M
Etrue = solvekep(Mtrue)

#setting up rotation matrices
normalpos = np.array([a*cos(Etrue)-a*e, a*sqrt(1-e**2)*sin(Etrue), 0])
rotateO = np.array([[cos(Oprime), -sin(Oprime), 0], [sin(Oprime), cos(Oprime), 0], [0,0,1]])
rotatei = np.array([[1,0,0],[0,cos(iprime), -sin(iprime)], [0,sin(iprime), cos(iprime)]])
rotatew = np.array([[cos(wprime), -sin(wprime), 0], [sin(wprime), cos(wprime), 0], [0,0,1]])

#calculate the ecliptic coordinates for the asteroid
r1ecliptic = rotateO @ rotatei @ rotatew @ normalpos

#put asteroid into vpython
asteroid = sphere(pos=vector(r1ecliptic[0], r1ecliptic[1], r1ecliptic[2])*150, radius=(15), color=color.white)
asteroid.trail = curve(color=color.white)
sun = sphere(pos=vector(0,0,0), radius=(50), color=color.yellow)

while (True):
    if time < 10:
        print(r1ecliptic)
    rate(150)
    time = time + 1

    #readjust M and E based on new time
    Mtrue = 2*pi/period*(time) + M
    Etrue = solvekep(Mtrue)

    #recalculate rectangular pos
    normalpos = np.array([a*cos(Etrue)-a*e, a*sqrt(1-e**2)*sin(Etrue), 0])

    #recalculate ecliptic coordinates
    r1ecliptic = rotateO @ rotatei @ rotatew @ normalpos

    #update asteroid position
    asteroid.pos = vector(r1ecliptic[0], r1ecliptic[1], r1ecliptic[2])*150
    asteroid.trail.append(pos=asteroid.pos)  
