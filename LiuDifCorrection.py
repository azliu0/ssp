from odlib import *
from MoGlib import *
from math import *
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
from astropy.io import fits


# these matrices are used to convert back and forth between ecliptic and equatorial coordinates
rotate_eps_cw = np.array([[1,0,0], [0, cos(eps), sin(eps)], [0,-sin(eps), cos(eps)]])
rotate_eps_ccw = np.array([[1,0,0], [0, cos(eps), -sin(eps)], [0, sin(eps), cos(eps)]])

ts, ras, decs, Rs = [], [], [], []
obs_file = open("2021MoGtestinput.txt", "r")

for line in obs_file:
    obs = line.split()
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
ra_1_og, ra_2_og, ra_3_og = radians(ras[0]), radians(ras[1]), radians(ras[2])
dec_1_og, dec_2_og, dec_3_og = radians(decs[0]), radians(decs[1]), radians(decs[2])
t_1_og, t_2_og, t_3_og = juliandate(ts[0][0], ts[0][1], ts[0][2], ts[0][3]), juliandate(ts[1][0], ts[1][1], ts[1][2], ts[1][3]), juliandate(ts[2][0], ts[2][1], ts[2][2], ts[2][3])


# R_1, R_2, R_3 = rotate_eps_ccw @ R_1, rotate_eps_ccw @ R_2, rotate_eps_ccw @ R_3

'''
--------------------------------------------------------------------------------------------------------------------------------
this section of the code represents the main method of gauss function
it takes the givens from the previous section to calculate r2 and r2hat using iteration from the method of gauss
--------------------------------------------------------------------------------------------------------------------------------
'''

# setting up known constants

rhohat_1,rhohat_2,rhohat_3 = sph_to_cart(ra_1_og, dec_1_og), sph_to_cart(ra_2_og, dec_2_og), sph_to_cart(ra_3_og, dec_3_og)

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
r_2_mag = r_2_mag_candidates[0]

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
while (abs(rho_2 - rho_2_prev) > 1e-12):
    rho_2_prev = rho_2

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

    t_1 = t_1 - rho_1/cAU
    t_2 = t_2 - rho_2/cAU
    t_3 = t_3 - rho_3/cAU

# converting to ecliptic
r_2_ec = rotate_eps_cw @ r_2
r_2_dot_ec = rotate_eps_cw @ r_2_dot

# orbital elements 
a,e,i, Omega, omega, M, E = orbital_elements_with_E(r_2_ec, r_2_dot_ec)


print(r_1, r_2, r_3)
print(np.linalg.norm(r_2))



# converting to ecliptic
r_2_ec = rotate_eps_cw @ r_2
r_2_dot_ec = rotate_eps_cw @ r_2_dot

def diff_ra_dec(r, R, delta):
    new_rho = r+R
    new_rho = new_rho / (np.linalg.norm(new_rho))

    r_x_og, r_y_og, r_z_og = new_rho[0], new_rho[1], new_rho[2]

    # modify x component 
    r_x, r_y, r_z = r_x_og - delta, r_y_og, r_z_og
    new_dec_1 = asin(r_z)
    new_ra_1 = quadrant_check(r_y/cos(new_dec_1), r_x/cos(new_dec_1))

    r_x, r_y, r_z = r_x_og + delta, r_y_og, r_z_og
    new_dec_2 = asin(r_z)
    new_ra_2 = quadrant_check(r_y/cos(new_dec_2), r_x/cos(new_dec_2))

    diff_ra_x = (new_ra_2 - new_ra_1)/(2*delta)
    diff_dec_x = (new_dec_2 - new_dec_1)/(2*delta)

    # modify y component 
    r_x, r_y, r_z = r_x_og, r_y_og - delta, r_z_og
    new_dec_1 = asin(r_z)
    new_ra_1 = quadrant_check(r_y/cos(new_dec_1), r_x/cos(new_dec_1))

    r_x, r_y, r_z = r_x_og, r_y_og + delta, r_z_og
    new_dec_2 = asin(r_z)
    new_ra_2 = quadrant_check(r_y/cos(new_dec_2), r_x/cos(new_dec_2))

    diff_ra_y = (new_ra_2 - new_ra_1)/(2*delta)
    diff_dec_y = (new_dec_2 - new_dec_1)/(2*delta)

    # modify z component 
    r_x, r_y, r_z = r_x_og, r_y_og, r_z_og - delta
    new_dec_1 = asin(r_z)
    new_ra_1 = quadrant_check(r_y/cos(new_dec_1), r_x/cos(new_dec_1))

    r_x, r_y, r_z = r_x_og, r_y_og, r_z_og + delta
    new_dec_2 = asin(r_z)
    new_ra_2 = quadrant_check(r_y/cos(new_dec_2), r_x/cos(new_dec_2))

    diff_ra_z = (new_ra_2 - new_ra_1)/(2*delta)
    diff_dec_z = (new_dec_2 - new_dec_1)/(2*delta)

    return diff_ra_x, diff_dec_x, diff_ra_y, diff_dec_y, diff_ra_z, diff_dec_z








