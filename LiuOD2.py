import numpy as np
from math import degrees
from odlib import k, retrieve_data, ang_momentum, orbital_elements

def main():
    # position and velocity vectors from ephemeris
    pos,velocity = retrieve_data("LiuInput.txt", "2018-Jul-14", "00:00:00.0000")

    # actual values for orbital elements
    t_a = 1.056800055745494E+00 
    t_e= 3.442331093323161E-01 
    t_i= 2.515525166662531E+01
    t_Omega= 2.362379806551942E+02 
    t_omega = 2.555046142766286E+02 
    t_nu_0= 1.589559248724372E+02
    t_period = 3.968146231808843E+02
    t_mean_anomaly= 1.404194576239256E+02 

    # calculated values for orbital elements
    a, e, i, Omega, omega, nu_0 = orbital_elements(pos, velocity)

    # errors
    er_a = 100*(abs(a - t_a)/t_a)
    er_e = 100*(abs(e - t_e)/t_e)
    er_i = 100*(abs(i - t_i)/t_i)
    er_Omega = 100*(abs(Omega - t_Omega)/t_Omega)
    er_omega = 100*(abs(omega - t_omega)/t_omega)
    er_nu_0 = 100*(abs(nu_0 - t_nu_0)/t_nu_0)

    #formatting print statements
    print("------------------------------------------------------------")
    print("semimajor axis:")
    print("expected =", t_a, "AU")
    print("calculated =", a, "AU")
    print("percent error =", er_a, "%")
    print("-----------------------------------------------")
    print("eccentricity:")
    print("expected =", t_e, "AU")
    print("calculated =", e, "AU")
    print("percent error =", er_e, "%")
    print("-----------------------------------------------")
    print("inclination:")
    print("expected =", t_i, "deg")
    print("calculated =", i, "deg")
    print("percent error =", er_i, "%")
    print("-----------------------------------------------")
    print("longitude of the ascending node:")
    print("expected =", t_Omega, "deg")
    print("calculated =", Omega, "deg")
    print("percent error =", er_Omega, "%")
    print("-----------------------------------------------")
    print("argument of perihelion:")
    print("expected =", t_omega, "deg")
    print("calculated =", omega, "deg")
    print("percent error =", er_omega, "%")
    print("------------------------------------------------------------")

main()