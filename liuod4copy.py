# Andrew Liu 
# OD4 Code
# SSP July 12th, 2021

from LSPR import dec_to_sg, ra_to_sg
import numpy as np
from math import pi, sqrt, cos, sin, asin, radians, degrees
from odlib import k, retrieve_data, quadrant_check, orbital_elements, juliandate, nr_kepler

# NOTE: code passes test cases perfectly

# given constants
sun = np.array([-3.881336047533506E-01, 8.619617590425438E-01, 3.736284118981542E-01])
eps = radians(23.4374)
# ra_check, dec_check = 15*(17+42/60+20.87/3600), 31+52/60+29/3600 #for testing

def ephemeris(pos, velocity, sun, t_2, t_1):
    a, e, i, Omega, omega, M_2 = orbital_elements(pos, velocity)

    print(a,e,i,Omega,omega,M_2)
    i, Omega, omega, M_2 = radians(i), radians(Omega), radians(omega), radians(M_2)
   

    # calculating elliptical eccentricity and mean anomaly at new date
    M = k*a**(-1.5)*(t_2 - t_1) + M_2
    print("M", M)
    E = nr_kepler(M, e)
    print("E", E)

    # rotating orbital plane position vector and sun vector into equatorial coordinates
    pos_op = np.array([a*cos(E)-a*e, a*sqrt(1-e**2)*sin(E), 0])

    print("pos_op", pos_op)
    # print("pos_op", pos_op)
    rotateO = np.array([[cos(Omega), -sin(Omega), 0], [sin(Omega), cos(Omega), 0], [0,0,1]])
    rotatei = np.array([[1,0,0],[0,cos(i), -sin(i)], [0,sin(i), cos(i)]])
    rotatew = np.array([[cos(omega), -sin(omega), 0], [sin(omega), cos(omega), 0], [0,0,1]])
    rotate_eps = np.array([[1,0,0], [0, cos(eps), -sin(eps)], [0,sin(eps), cos(eps)]])

    pos_eq = rotate_eps @ rotateO @ rotatei @ rotatew @ pos_op
    sun_eq = rotate_eps @ sun

    print("sun_eq", sun_eq)

    # print("asdfasdf", pos_eq)
    # print("sun_eq", sun_eq)

    # calculating RA and dec
    rho = pos_eq + sun_eq
    rho = rho / (np.linalg.norm(rho))

    print("rho", rho)

    dec = asin(rho[2])
    ra = quadrant_check(rho[1]/cos(dec), rho[0]/cos(dec))

    ra, dec = degrees(ra), degrees(dec)
    print(ra,dec)

    print("asdfasdf", ra, dec)
    print("asdfasdf", ra_to_sg(ra), dec_to_sg(dec))

sun_4 = np.array([-4.491811083018563E-01, 9.115858819358923E-01, -2.895047278592023E-05])
sun_2 = np.array([-2.621569363571260E-01, 9.823620200990218E-01, -8.140802287552636E-05])
r2 = np.array([0.37736089, -1.52293765, 0.07793864])
r2dot = np.array([0.78060228, 0.55026231, 0.02776518])
t_2, t_3, t_4 = 2459402.6245254627, 2459408.6705439813, 2459414.4610069445

sun_sanjana = np.array([-2.346938368583631E-01,9.892703177219488E-01,-4.710393708209949E-05])
t_2 = juliandate(2021, 7, 24, 0)
t_3 = juliandate(2021, 7, 5, 10.8)
r2 = np.array([4.034081027772663E-01, -1.284498396379768E+00, 5.641805193316963E-01])
r2dot = np.array( [0.8113729348954035, -0.02441906025580499, 0.2499992985301538])



ephemeris(r2, r2dot, sun_sanjana, t_2, t_3)



