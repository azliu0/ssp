from astropy.units.function.logarithmic import M_bol
from odlib import *
from MoGlib import *
from math import *
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
from astropy.io import fits

def uncertainty():
    # these matrices are used to convert back and forth between ecliptic and equatorial coordinates
    rotate_eps_cw = np.array([[1,0,0], [0, cos(eps), sin(eps)], [0,-sin(eps), cos(eps)]])
    rotate_eps_ccw = np.array([[1,0,0], [0, cos(eps), -sin(eps)], [0, sin(eps), cos(eps)]])

    ts, ras, decs, Rs = [], [], [], []
    obs_file = open("LiuODinput_actual.txt", "r")

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

    R_1, R_2, R_3 = rotate_eps_ccw @ R_1, rotate_eps_ccw @ R_2, rotate_eps_ccw @ R_3

    ra_sigma_2, dec_sigma_2 = RMS_uncertainty("corr_obs2.fits")
    ra_sigma_3, dec_sigma_3 = RMS_uncertainty("corr_obs3.fits")
    ra_sigma_4, dec_sigma_4 = RMS_uncertainty("corr_obs4.fits")

    num_reps = 0
    print("------------------------------------------------")
    print("how many iterations would you like to simulate (recommended >=1000)?")
    x = input()
    total_reps = int(x)
    a_s, e_s, i_s, omega_s, Omega_s, M_s = [], [], [], [], [], []

    print("------------------------------------------------")
    print("starting simulation...")

    while num_reps < total_reps:
        ra_1, ra_2, ra_3 = np.random.normal(ra_1_og,ra_sigma_2), np.random.normal(ra_2_og,ra_sigma_3), np.random.normal(ra_3_og,ra_sigma_4)
        dec_1, dec_2, dec_3 = np.random.normal(dec_1_og,dec_sigma_2), np.random.normal(dec_2_og,dec_sigma_3), np.random.normal(dec_3_og,dec_sigma_4)

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
        r_2_mag = r_2_mag_candidates[0]
        rho_2_mag = rho_2_candidates[0]

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
            count += 1

        # converting to ecliptic
        r_2_ec = rotate_eps_cw @ r_2
        r_2_dot_ec = rotate_eps_cw @ r_2_dot

        # orbital elements 
        a,e,i, Omega, omega, M, E = orbital_elements_with_E(r_2_ec, r_2_dot_ec)

        a_s.append(a)
        e_s.append(e)
        i_s.append(i)
        Omega_s.append(Omega)
        omega_s.append(omega)
        M_s.append(M)

        num_reps += 1
        if (num_reps % 500 == 0):
            print("iterations complete:", num_reps)

    print("------------------------------------------------")
    input("Simulation complete! Press enter to display results...")

    num_bins = 101

    def hist(title, x_label, units, el, jpl_el):
        N, bins, patches = plt.hist(el, num_bins, facecolor='pink')
        #  edgecolor = 'black'
        bin_centers = (bins[:-1] + bins[1:])/2
        col = abs(bin_centers - np.mean(el))
        fracs = 1 - col / max(col)
        for frac, patch in zip(fracs, patches):
            color = plt.cm.RdPu(frac)
            patch.set_facecolor(color)
        plt.axvline(x=np.mean(el), color = 'lightsalmon', label = r'$\mu = $'+str(round(np.mean(el),3)) + " " + units)
        plt.axvline(x=jpl_el, color = 'red', label = r'$\mu$'+' (JPL) = '+str(round(jpl_el,3)) + " " + units)
        plt.axvline(x=np.mean(el)+np.std(el), color = 'lightsalmon', linestyle = 'dotted', label = r'$\sigma =$'+str(round(np.std(el),5)) + " " + units)
        plt.axvline(x=np.mean(el)-np.std(el), color = 'lightsalmon', linestyle = 'dotted')
        plt.axvline(x=np.mean(el)+2*np.std(el), color = 'lightsalmon', linestyle = 'dotted')
        plt.axvline(x=np.mean(el)-2*np.std(el), color = 'lightsalmon', linestyle = 'dotted')

        plt.xlabel(x_label)
        plt.ylabel('# of instances')
        plt.title(title)
        plt.legend(loc="upper right")

    # jpl_a = 2.726573219901208
    # jpl_e = 0.531749369668808
    # jpl_i = 3.8752
    # jpl_Omega = 238.836
    # jpl_omega = 107.9214
    # jpl_M = 339.74

    # EC= 5.363327337998935E-01 QR= 1.259944760903456E+00 IN= 3.887497846537230E+00
    #  OM= 2.383885689081052E+02 W = 1.083143069347813E+02 Tp=  2459492.572394193150
    #  N = 2.200323168118316E-01 MA= 3.415388815130426E+02 TA= 2.974294498648268E+02
    #  A = 2.717346797476311E+00 AD= 4.174748834049167E+00 PR= 1.636123298687378E+03

    # EC= 5.363327337998935E-01 QR= 1.259944760903456E+00 IN= 3.887497846537230E+00
    #  OM= 2.383885689081052E+02 W = 1.083143069347813E+02 Tp=  2459492.572394193150
    #  N = 2.200323168118316E-01 MA= 3.415388815130426E+02 TA= 2.974294498648268E+02
    #  A = 2.717346797476311E+00 AD= 4.174748834049167E+00 PR= 1.636123298687378E+03

    jpl_a = 2.717346797476311E+00 
    jpl_e = 5.363327337998935E-01
    jpl_i = 3.887497846537230E+00
    jpl_Omega = 2.383885689081052E+02
    jpl_omega = 1.083143069347813E+02
    jpl_M = 3.415388815130426E+02

    mean_a, mean_e, mean_i, mean_Omega, mean_omega, mean_M = np.mean(a_s), np.mean(e_s), np.mean(i_s), np.mean(Omega_s), np.mean(omega_s), np.mean(M_s)
    std_a, std_e, std_i, std_Omega, std_omega, std_M = np.std(a_s), np.std(e_s), np.std(i_s), np.std(Omega_s), np.std(omega_s), np.std(M_s)

    print("------------------------------------------------------------")
    print("semimajor axis:")
    print("expected =", jpl_a, "AU")
    print("calculated mean value (+/- SD)=", str(mean_a) + "+/-" + str(std_a), "AU")
    print("percent error =", error(jpl_a, mean_a), "%")
    print("-----------------------------------------------")
    print("eccentricity:")
    print("expected =", jpl_e, "AU")
    print("calculated mean value (+/- SD)=", str(mean_e) + "+/-" + str(std_e))
    print("percent error =", error(jpl_e, mean_e), "%")
    print("-----------------------------------------------")
    print("inclination:")
    print("expected =", jpl_i, "deg")
    print("calculated mean value (+/- SD)=", str(mean_i) + "+/-" + str(std_i), "deg")
    print("percent error =", error(jpl_i, mean_i), "%")
    print("-----------------------------------------------")
    print("longitude of the ascending node:")
    print("expected =", jpl_Omega, "deg")
    print("calculated mean value (+/- SD)=", str(mean_Omega) + "+/-" + str(std_Omega), "deg")
    print("percent error =", error(jpl_Omega, mean_Omega), "%")
    print("-----------------------------------------------")
    print("argument of perihelion:")
    print("expected =", jpl_omega, "deg")
    print("calculated mean value (+/- SD)=", str(mean_omega) + "+/-" + str(std_omega), "deg")
    print("percent error =", error(jpl_omega, mean_omega), "%")
    print("-----------------------------------------------")
    print("mean anomaly (M):")
    print("expected =", jpl_M, "deg")
    print("calculated mean value (+/- SD)=", str(mean_M) + "+/-" + str(std_M), "deg")
    print("percent error =", error(jpl_M, mean_M), "%")
    print("-----------------------------------------------")

    plot_a = plt.figure(1)
    hist(str(total_reps) + " iterations of Semi-Major Axis (a)", "Semi-Major Axis (AU)", "AU", a_s, jpl_a)

    plot_e = plt.figure(2)
    hist(str(total_reps) + " iterations iterations of eccentricity (e)", "Eccentricity", "", e_s, jpl_e)

    plot_i = plt.figure(3)
    hist(str(total_reps) + " iterations of inclination (i)", "Inclination (degrees)", "degrees", i_s, jpl_i)

    plot_O = plt.figure(4)
    hist(str(total_reps) + ' iterations of Longitude of Ascending Node ('+r'$\Omega$' + ')', "Longitude of Ascending Node (degrees)", "degrees", Omega_s, jpl_Omega)

    plot_w = plt.figure(5)
    hist(str(total_reps) + ' iterations of Argument of Perihelion ('+r'$\omega$' + ')', "Argument of Perihelion (degrees)", "degrees", omega_s, jpl_omega)

    plot_M = plt.figure(6)
    hist(str(total_reps) + ' iterations of Mean Anomaly (M)', "Mean Anomaly (degrees)", "degrees", M_s, jpl_M)

    print("------------------------------------------------")
    print("To continue to run this program, close all histogram windows:")

    plt.show()







