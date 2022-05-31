
import numpy as np
from math import *
from odlib import *


def SEL(taus, sun2, rhohat2, Ds):
    roots, rhos = [], []

    tau_1,tau_3, tau = taus[0],taus[1],taus[2]
    D_0, D_21, D_22, D_23 = Ds[0], Ds[1], Ds[2], Ds[3]
    
    A_1 = tau_3/tau
    A_3 = -tau_1/tau
    B_1 = A_1/6 * (tau**2 - tau_3**2)
    B_3 = A_3/6*(tau**2 - tau_1**2)

    A = (A_1*D_21 - D_22 + A_3*D_23) / (-D_0)
    B = (B_1*D_21 + B_3*D_23)/ (-D_0)
    E = -2*(np.dot(rhohat2, sun2))
    F = np.linalg.norm(sun2)**2

    a = -(A**2 + A*E + F)
    b = -(2*A*B + B*E)
    c = -B**2

    for num in np.polynomial.polynomial.polyroots([c, 0, 0, b, 0, 0, a, 0, 1]):
        if abs(num.imag) < 1e-5 and num.real > 0:
            roots.append(num.real)
            rhos.append(A + B/(num.real**3))
    
    return roots, rhos

def f0g0(tau, r2_mag):
    u = 1/r2_mag**3
    f = 1 - 1/2*u*tau**2
    g = tau - 1/6*u*tau**3
    return f,g

def iterate_deltaE(tau, r2, r2dot):
    r2_mag, r2dot_mag = np.linalg.norm(r2), np.linalg.norm(r2dot)
    a = 1 / (2/r2_mag - r2dot_mag**2)
    tol = 1e-12
    n = sqrt(1/a**3)
    delta_E, delta_E_prev = n*tau, -1
    def func(x):
        return x - (1 - r2_mag/a)*sin(x)+(np.dot(r2, r2dot))/(n*a**2)*(1-cos(x))-n*tau
    def func_d(x):
        return 1 - (1 - r2_mag/a)*cos(x)+(np.dot(r2, r2dot))/(n*a**2)*(sin(x))
    while (abs(delta_E - delta_E_prev) > tol):
        delta_E_prev = delta_E
        delta_E = delta_E - func(delta_E)/func_d(delta_E)
    f,g = 1 - a/r2_mag * (1-cos(delta_E)), tau + 1/n*(sin(delta_E)-delta_E)
    return f,g

def iterate_taylor(tau, r2, r2dot):
    r2_mag = np.linalg.norm(r2)
    u = 1/r2_mag**3
    z = np.dot(r2, r2dot)/(r2_mag**2)
    q = np.dot(r2dot, r2dot)/(r2_mag**2) - u
    f = 1 - 1/2*u*tau**2 + 1/2*u*z*tau**3 + 1/24*(3*u*q-15*u*z**2+u**2)*tau**4
    g = tau - 1/6*u*tau**3 + 1/4*u*z*tau**4
    return f,g