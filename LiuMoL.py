from odlib import *
from MoGlib import *
from math import *
import numpy as np
from LiuODvisual import *

rotate_eps_cw = np.array([[1,0,0], [0, cos(eps), sin(eps)], [0,-sin(eps), cos(eps)]])
rotate_eps_ccw = np.array([[1,0,0], [0, cos(eps), -sin(eps)], [0, sin(eps), cos(eps)]])

# reading in given data from input file

# ts, ras, decs, Rs = [], [], [], []
# obs_file = open("chuiinput.txt", "r")

# for line in (obs_file):
#     obs = line.split()
#     for i in range(len(obs)):
#         if i == 3:
#             temp = obs[i].split(":")
#             obs[i] = float(temp[0]) + float(temp[1])/60 + float(temp[2])/3600
#         if i == 4:
#             obs[i] = ra_to_degrees(obs[i])
#         if i == 5:
#             obs[i] = dec_to_degrees(obs[i])
#         else:
#             obs[i] = float(obs[i])
    
#     ts.append(obs[:4])
#     ras.append(obs[4])
#     decs.append(obs[5])
#     Rs.append(obs[6:])

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
print("which observations would you like to use (3 integers separated by commas):")
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

print("decs", decs)
print("ras", ras)
R_1, R_2, R_3 = np.array(Rs[0]), np.array(Rs[1]), np.array(Rs[2]) 
ra_1, ra_2, ra_3 = radians(ras[0]), radians(ras[1]), radians(ras[2])
dec_1, dec_2, dec_3 = radians(decs[0]), radians(decs[1]), radians(decs[2])
t_1_og, t_2_og, t_3_og = juliandate(ts[0][0], ts[0][1], ts[0][2], ts[0][3]), juliandate(ts[1][0], ts[1][1], ts[1][2], ts[1][3]), juliandate(ts[2][0], ts[2][1], ts[2][2], ts[2][3])

R_1, R_2, R_3 = rotate_eps_ccw @ R_1, rotate_eps_ccw @ R_2, rotate_eps_ccw @ R_3

tau_1, tau_3 = k*(t_1_og-t_2_og), k*(t_3_og-t_2_og)
tau = tau_3-tau_1

rhohat_1,rhohat_2,rhohat_3 = sph_to_cart(ra_1, dec_1), sph_to_cart(ra_2, dec_2), sph_to_cart(ra_3, dec_3)

rhohat_2_dot = (tau_3**2*(rhohat_1 - rhohat_2) - tau_1**2 * (rhohat_3-rhohat_2))/(tau_1*tau_3*tau)
rhohat_2_ddot = -2*(tau_3 * (rhohat_1 - rhohat_2) - tau_1*(rhohat_3 - rhohat_2))/(tau_1*tau_3*tau)

print(rhohat_2_dot, rhohat_2_ddot)

A = np.dot(np.cross(rhohat_2, rhohat_2_dot), R_2)/(np.dot(np.cross(rhohat_2, rhohat_2_dot), rhohat_2_ddot))
B = (1 + 1/328900.5)*A/(np.linalg.norm(R_2)**3)


D_0 = np.dot(rhohat_1, np.cross(rhohat_2, rhohat_3))
D_21 = np.dot(np.cross(rhohat_1, R_1), rhohat_3)
D_22 = np.dot(np.cross(rhohat_1, R_2), rhohat_3)
D_23 = np.dot(np.cross(rhohat_1, R_3), rhohat_3)

taus, Ds = [tau_1, tau_3, tau], [D_0, D_21, D_22, D_23]

print("taus", taus)
print("Ds", Ds)
r_2_mag_candidates, rho_2_candidates = SEL(taus, R_2, rhohat_2, Ds)

for i in range(len(r_2_mag_candidates)-1,-1,-1):
    if r_2_mag_candidates[i] < 0 or rho_2_candidates[i] < 0:
        r_2_mag_candidates.pop(i)
        rho_2_candidates.pop(i)

print("\ncandidates for r_2:")
print("------------------------------------------------------------------------------------------")
for i in range(len(r_2_mag_candidates)):
    print("option", str(i+1) + ":", "r_2 =", str(r_2_mag_candidates[i]), "AU; corresponding rho_2 =", str(rho_2_candidates[i]), "AU")
print("------------------------------------------------------------------------------------------\n")

def is_int(x):
        try:
            x = int(x)
        except ValueError:
            return False
        return True

print("choose which option number:")
x = input()
while (not is_int(x) or int(x) not in range(1, len(r_2_mag_candidates)+1)):
    print("invalid input, please choose which option number:")
    x = input()
x = int(x)
r_2_mag = r_2_mag_candidates[x-1]
rho_2_mag = rho_2_candidates[x-1]
# rho_2_mag = A / r_2_mag**3 - B

print(r_2_mag, rho_2_mag)
print(A,B)

r_2_mag, rho_2_mag = 1.0, 0.0

r_2_prev, rho_2_prev = -1,-1

count = 0
while (abs(r_2_prev - r_2_mag) > 1e-12 or abs(rho_2_prev - rho_2_mag) > 1e-12):
    r_2_prev, rho_2_prev = r_2_mag, rho_2_mag
    r_2_mag = sqrt(rho_2_prev**2 + np.linalg.norm(R_2)**2 - 2*np.dot(R_2, rho_2_prev*rhohat_2))
    rho_2_mag = A / r_2_mag**3 - B
    print(count, r_2_mag, rho_2_mag)
    count += 1

print(r_2_mag, rho_2_mag)

rho_2_mag_dot = -1/2*(1/r_2_mag**3 - (1 + 1/328900.56)/(np.linalg.norm(R_2)**3))*(np.dot(np.cross(rhohat_2, rhohat_2_ddot), R_2))/(np.dot(np.cross(rhohat_2, rhohat_2_dot), rhohat_2_ddot))

R_2_dot = (R_3 - R_1)/tau

r_2 = rho_2_mag*rhohat_2 - R_2
r_2_dot = rho_2_mag_dot*rhohat_2 + rho_2_mag*rhohat_2_dot - R_2_dot

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
# print("\n with rho2 =", rho_2, "AU")
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