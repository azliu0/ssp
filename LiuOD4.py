# Andrew Liu 
# OD4 Code
# SSP July 12th, 2021

import numpy as np
from math import pi, sqrt, cos, sin, asin, radians, degrees
from odlib import k, retrieve_data, quadrant_check, orbital_elements, juliandate, nr_kepler

# NOTE: code passes test cases perfectly

# given constants
pos,velocity = retrieve_data("LiuInput_OD4.txt", "2018-Jul-14", "00:00:00.0000")
sun = np.array([-6.573717007738460E-01, 7.730214598589336E-01, -6.405807721620831E-05])
eps = radians(23.4374)
ra_check, dec_check = 15*(17+42/60+20.87/3600), 31+52/60+29/3600 #for testing

def main():
    a, e, i, Omega, omega, M_2 = orbital_elements(pos, velocity)
    print(a,e, i, Omega, omega, M_2)
    i, Omega, omega, M_2 = radians(i), radians(Omega), radians(omega), radians(M_2)
    print(a,e, i, Omega, omega, M_2)
   

    # calculating elliptical eccentricity and mean anomaly at new date
    M = k*a**(-1.5)*(juliandate(2018, 8, 3, 0) - juliandate(2018, 7, 14, 0)) + M_2
    print(k*a**(-1.5)*(juliandate(2018, 8, 3, 0) - juliandate(2018, 7, 14, 0)))

    print(k*a**(-1.5))
    print(0.01720209895**2)
    print(M)
    print(degrees(M))
    E = nr_kepler(M, e)
    print("E,", E)

    # rotating orbital plane position vector and sun vector into equatorial coordinates
    pos_op = np.array([a*cos(E)-a*e, a*sqrt(1-e**2)*sin(E), 0])
    print("pos_op", pos_op)
    rotateO = np.array([[cos(Omega), -sin(Omega), 0], [sin(Omega), cos(Omega), 0], [0,0,1]])
    rotatei = np.array([[1,0,0],[0,cos(i), -sin(i)], [0,sin(i), cos(i)]])
    rotatew = np.array([[cos(omega), -sin(omega), 0], [sin(omega), cos(omega), 0], [0,0,1]])
    rotate_eps = np.array([[1,0,0], [0, cos(eps), -sin(eps)], [0,sin(eps), cos(eps)]])

    pos_eq = rotate_eps @ rotateO @ rotatei @ rotatew @ pos_op
    sun_eq = rotate_eps @ sun

    print("asdfasdf", pos_eq)
    print("sun_eq", sun_eq)

    # calculating RA and dec
    rho = pos_eq + sun_eq
    rho = rho / (np.linalg.norm(rho))
    dec = asin(rho[2])
    ra = quadrant_check(rho[1]/cos(dec), rho[0]/cos(dec))


    print("-----------------------------------------------")
    print("RA:")
    print("expected =", ra_check, "deg")
    print("calculated =", degrees(ra), "deg")
    print("percent error =", abs(ra_check - degrees(ra))/ra_check * 100, "%")
    print("-----------------------------------------------")
    print("Dec:")
    print("expected =", dec_check, "deg")
    print("calculated =", degrees(dec), "deg")
    print("percent error =", abs(dec_check - degrees(dec))/dec_check * 100, "%")
    print("-----------------------------------------------")

main()

