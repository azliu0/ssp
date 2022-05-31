from odlib import *
from MoGlib import *
from math import *
import numpy as np
from LiuODvisual import *
from LiuODUncertainty import *
import sys

# NOTE: code works exactly right

# these matrices are used to convert back and forth between ecliptic and equatorial coordinates
rotate_eps_cw = np.array([[1,0,0], [0, cos(eps), sin(eps)], [0,-sin(eps), cos(eps)]])
rotate_eps_ccw = np.array([[1,0,0], [0, cos(eps), -sin(eps)], [0, sin(eps), cos(eps)]])

'''
--------------------------------------------------------------------------------------------------------------------------------
this section of the code reads in the data from the input file. 
from the input file we are given time of observations, RAs, DECs, and corresponding E-S vectors 
the code takes in values and converts them into the values useful for method of gauss: rhohats and taus
--------------------------------------------------------------------------------------------------------------------------------
'''

ts, ras, decs, Rs = [], [], [], []
obs_file = open("LiuODinput.txt", "r")
obs_data = []

for line in (obs_file):
    obs_data += [line.split()]

# asking user input to choose which observations they would like to use
print("----------------------------------------------------------------")
print("you have", len(obs_data), "observations available from the following dates:\n")
for i in range(len(obs_data)):
    print("observation", str(i+1) + ":", obs_data[i][0] + "/" + obs_data[i][1] + "/" + obs_data[i][2], "at", obs_data[i][3], "UTC")
print("----------------------------------------------------------------")
print("which observations would you like to use? (3 integers separated by commas):")
x = input()
x = x.split(",")

# error checking user entry, making sure that they are three valid observations
while (len(x) != 3 or not is_int(x[0]) or not is_int(x[1]) or not is_int(x[2]) or int(x[0]) not in range(len(obs_data)+1) or int(x[1]) not in range(len(obs_data)+1) or int(x[2]) not in range(len(obs_data)+1)):
    if (len(x) != 3 or not is_int(x[0]) or not is_int(x[1]) or not is_int(x[2])):
        print("----------------------------------------------------------------")
        print("(ERROR) invalid input, please enter 3 integers separated by commas:")
    else:
        print("----------------------------------------------------------------")
        print("(ERROR) you only have the following observations available:")
        for i in range(len(obs_data)):
            print(i+1, end = " ")
        print()
        print("----------------------------------------------------------------")
        print("please rechoose three from the list above only:")  
    x = input()
    x = x.split(",")

for i in range(len(x)):
    x[i] = int(x[i])-1
x.sort()

for i in range(len(obs_data)-1, -1, -1):
    if i not in x:
        obs_data.pop(i)

# finally, reading in data from the three desired observations
for obs in obs_data:
    for i in range(len(obs)):
        if i == 3:
            temp = obs[i].split(":")
            obs[i] = float(temp[0]) + float(temp[1])/60 + float(temp[2])/3600
        if i == 4:
            obs[i] = ra_to_degrees(obs[i])
        if i == 5:
            obs[i] = dec_to_degrees(obs[i])
        else:
            obs[i] = float(obs[i])

    ts.append(obs[:4])
    ras.append(obs[4])
    decs.append(obs[5])
    Rs.append(obs[6:])

R_1, R_2, R_3 = np.array(Rs[0]), np.array(Rs[1]), np.array(Rs[2]) 
ra_1, ra_2, ra_3 = radians(ras[0]), radians(ras[1]), radians(ras[2])
dec_1, dec_2, dec_3 = radians(decs[0]), radians(decs[1]), radians(decs[2])
t_1_og, t_2_og, t_3_og = juliandate(ts[0][0], ts[0][1], ts[0][2], ts[0][3]), juliandate(ts[1][0], ts[1][1], ts[1][2], ts[1][3]), juliandate(ts[2][0], ts[2][1], ts[2][2], ts[2][3])

R_1, R_2, R_3 = rotate_eps_ccw @ R_1, rotate_eps_ccw @ R_2, rotate_eps_ccw @ R_3
'''
--------------------------------------------------------------------------------------------------------------------------------
this section of the code represents the main method of gauss function
it takes the givens from the previous section to calculate r2 and r2hat using iteration from the method of gauss
--------------------------------------------------------------------------------------------------------------------------------
'''

# setting up known constants

rhohat_1,rhohat_2,rhohat_3 = sph_to_cart(ra_1, dec_1), sph_to_cart(ra_2, dec_2), sph_to_cart(ra_3, dec_3)

D_0 = np.dot(rhohat_1, np.cross(rhohat_2, rhohat_3))

D_11 = np.dot(np.cross(R_1, rhohat_2), rhohat_3)
D_21 = np.dot(np.cross(rhohat_1, R_1), rhohat_3)
D_31 = np.dot(rhohat_1, np.cross(rhohat_2, R_1))

D_12 = np.dot(np.cross(R_2, rhohat_2), rhohat_3)
D_22 = np.dot(np.cross(rhohat_1, R_2), rhohat_3)
D_32 = np.dot(rhohat_1, np.cross(rhohat_2, R_2))

D_13 = np.dot(np.cross(R_3, rhohat_2), rhohat_3)
D_23 = np.dot(np.cross(rhohat_1, R_3), rhohat_3)
D_33 = np.dot(rhohat_1, np.cross(rhohat_2, R_3))

# setting up intial values for iteration using the scalar equations of lagrange

tau_1, tau_3 = k*(t_1_og-t_2_og), k*(t_3_og-t_2_og)
tau = tau_3-tau_1

taus, Ds = [tau_1, tau_3, tau], [D_0, D_21, D_22, D_23]

r_2_mag_candidates, rho_2_candidates = SEL(taus, R_2, rhohat_2, Ds)

for i in range(len(r_2_mag_candidates)-1,-1,-1):
    if r_2_mag_candidates[i] < 0 or rho_2_candidates[i] < 0:
        r_2_mag_candidates.pop(i)
        rho_2_candidates.pop(i)

if (len(r_2_mag_candidates) == 0):
    print("----------------------------------------------------------------")
    print("no initial values available, please choose different observations. sorry ):")
    sys.exit()

# asking user entry for which roots to proceed with the main iteration
print("----------------------------------------------------------------")
print("You have the following candidates available for r_2:\n")
for i in range(len(r_2_mag_candidates)):
    print("option", str(i+1) + ":", "r_2 =", str(r_2_mag_candidates[i]), "AU; corresponding rho_2 =", str(rho_2_candidates[i]), "AU")
print("----------------------------------------------------------------")
print("please choose which candidate you would like to select (integer):")
x = input()
while (not is_int(x) or int(x) not in range(1, len(r_2_mag_candidates)+1)):
    print("invalid input, please choose a valid integer option number:")
    x = input()
x = int(x)
r_2_mag = r_2_mag_candidates[x-1]
rho_2_mag = rho_2_candidates[x-1]

# more initialization before main iteration loop

f_1, g_1 = f0g0(tau_1, r_2_mag)
f_3, g_3 = f0g0(tau_3, r_2_mag)

c_1, c_2, c_3, d_1, d_3 = g_3 / (f_1*g_3 - g_1*f_3), -1, -g_1 / (f_1*g_3 - g_1*f_3), -f_3 / (f_1*g_3 - g_1*f_3), f_1 / (f_1*g_3 - g_1*f_3), 
rho_1 = (c_1*D_11 + c_2*D_12+c_3*D_13)/(c_1*D_0)
rho_2 = (c_1*D_21 + c_2*D_22+c_3*D_23)/(c_2*D_0)
rho_3 = (c_1*D_31 + c_2*D_32+c_3*D_33)/(c_3*D_0)

r_1 = rhohat_1*rho_1 - R_1
r_2 = rhohat_2*rho_2 - R_2
r_3 = rhohat_3*rho_3 - R_3

r_2_dot = d_1*r_1 + d_3*r_3

t_1 = t_1_og - rho_1/cAU
t_2 = t_2_og - rho_2/cAU
t_3 = t_3_og - rho_3/cAU

tau_1, tau_3 = k*(t_1-t_2), k*(t_3-t_2)

# main iteration loop 

rho_2_prev = -1
count = 0
print("------------------")
input("Press enter to begin main iteration loop ...")

print("------------------")
print("\nMain Iteration Loop")
while (abs(rho_2 - rho_2_prev) > 1e-12):
    rho_2_prev = rho_2
    # tau_1, tau_3 = k*(t_1-t_2), k*(t_3-t_2)

    f_1, g_1 = iterate_deltaE(tau_1, r_2, r_2_dot)
    f_3, g_3 = iterate_deltaE(tau_3, r_2, r_2_dot)

    c_1, c_2, c_3, d_1, d_3 = g_3 / (f_1*g_3 - g_1*f_3), -1, -g_1 / (f_1*g_3 - g_1*f_3), -f_3 / (f_1*g_3 - g_1*f_3), f_1 / (f_1*g_3 - g_1*f_3), 
    
    rho_1 = (c_1*D_11 + c_2*D_12+c_3*D_13)/(c_1*D_0)
    rho_2 = (c_1*D_21 + c_2*D_22+c_3*D_23)/(c_2*D_0)
    rho_3 = (c_1*D_31 + c_2*D_32+c_3*D_33)/(c_3*D_0)

    r_1 = rhohat_1*rho_1 - R_1
    r_2 = rhohat_2*rho_2 - R_2
    r_3 = rhohat_3*rho_3 - R_3

    r_2_dot = d_1*r_1 + d_3*r_3

    t_1 = t_1_og - rho_1/cAU
    t_2 = t_2_og - rho_2/cAU
    t_3 = t_3_og - rho_3/cAU

    print(str(count+1).zfill(2) + ":", "change in rho2 =", rho_2 - rho_2_prev, "AU;", "light-travel time =", rho_2/cAU * 86400, "s.")
    count += 1

'''
--------------------------------------------------------------------------------------------------------------------------------
this section of the code utilizes the orbit determination code to print out the orbital elements
based off of the iterated r2 and r2hat values from the previous section.
this part of the code is mostly print statement formatting
--------------------------------------------------------------------------------------------------------------------------------
'''

print("------------------")
input("Iteration complete! Press enter to display results ...")
# converting to ecliptic
r_2_ec = rotate_eps_cw @ r_2
r_2_dot_ec = rotate_eps_cw @ r_2_dot

# orbital elements 
a,e,i, Omega, omega, M, E = orbital_elements_with_E(r_2_ec, r_2_dot_ec)
n = sqrt(1/a**3) #gauss units, radians
M_t = degrees(k*n*(juliandate(2021, 7, 24, 7) - t_2_og) + radians(M))

# print statements
print("\nIn", count, "iterations, r2 and r2dot converged to")
print("r2 =", r_2, "=", np.linalg.norm(r_2), "AU")
print("r2dot =", r_2_dot, "=", np.linalg.norm(r_2_dot)*k, "AU/day")
print("in cartesian equatorial coordinates.")
print("or")
print("r2 =", r_2_ec, "=", np.linalg.norm(r_2_ec), "AU")
print("r2dot =", r_2_dot_ec, "=", np.linalg.norm(r_2_dot_ec)*k, "AU/day")
print("in cartesian ecliptic coordinates.")
print("\n with rho2 =", rho_2, "AU")
print("\n ORBITAL ELEMENTS")
print("\t a =", a, "au")
print("\t e =", e)
print("\t i =", i, "degrees")
print("\t omega =", omega, "degrees")
print("\t Omega =", Omega, "degrees")
print("\t M =", M, "degrees at central obs. (JD =", str(t_2_og) + ")")
print("\t M =", M_t, "degrees at 07/24/2021 07:00 (JD =", str(juliandate(2021, 7, 24, 7)) + ")")
print("Also...")
print("\t E =", E, "degrees at central obs.")
print("\t n =", degrees(n)*k, "deg/day")
print("\t JD of last perihelion passage =", t_2_og - radians(M)/n*1/k)
print("\t P =", 2*pi/n * 1/k * 1/365.2425, "yrs")
print("\t P =", 2*pi/n * 1/k, "days")

print("------------------")
input("press enter to continue:")

print("------------------")
print("would you like to run a simulation to calculate uncertainties and SDs for each orbital element?\n This will also compare computed data vs. JPL data (y/n):")
x = input()
valid = ["yes", "y", "no", "n"]
while (x not in valid):
    print("invalid input. please enter (y/n) if you would like to run the simulation:")
    x = input()
while (x == "yes" or x == "y"):
    uncertainty()
    print("------------------")
    print("would you like to run another simulation? (y/n)")
    x = input()
    while (x not in valid):
        print("invalid input. please enter (y/n) if you would like to run the simulation:")
        x = input()

print("------------------")
input("press enter to continue:")

print("------------------")
print("Would you like to visualize the orbit of your asteroid? (y/n)")
x = input()
valid = ["yes", "y", "no", "n"]
while (x not in valid):
    print("invalid input. please enter (y/n) if you would like to visualize the orbit of your asteroid:")
    x = input()
if (x == "yes" or x == "y"):
    visualize()

print("------------------")
print("  ^~^  ,")
print(" ('Y') )")
print(" /   \/ ")
print("(\|||/) ")
print("------------------")
print("~thanks for running!")






