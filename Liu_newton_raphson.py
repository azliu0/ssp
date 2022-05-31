# Andrew Liu 
# SSP July 12th, 2021

# NOTE: code runs perfectly as expected

import numpy as np
from math import sin, cos
from scipy import optimize

M,e = 0.42, 0.8

#error function
def f(x):
    return x - e*sin(x)-M

#newton-raphson
def iterate_nr(E_guess):
    E_guess = E_guess - (M - (E_guess - e*sin(E_guess)))/(e*cos(E_guess) - 1)
    return E_guess

#scipy
def iterate_sc(E_guess):
    return optimize.newton(f, E_guess)

def main():
    M, e = 0.42, 0.8
    counter, tol = 0, 0.5e-5

    E_guess = M
    E_sc = iterate_sc(E_guess)
    E_guess_prev = -1
    while (abs(E_guess - E_guess_prev) > tol):
        E_guess_prev = E_guess
        E_guess = iterate_nr(E_guess)
        counter += 1
    
    print("E_homemade =", E_guess)
    print("E_scipy =", E_sc)
    print("Convergence parameter =", tol)
    print("Number of iterations: ", counter)

main()

