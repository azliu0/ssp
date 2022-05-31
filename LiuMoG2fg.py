# Andrew Liu MoG2: f&g
# SSP July 12th, 2021

# NOTE: code works exactly as expected

import numpy as np
from math import sin, cos

# flag = "3rd" if wanting 3rd order approximation 
#      = "4th" if wanting 4th order approximation
#      = tolerance if wanting an approximation using the exact formula
def fg(tau, r2_v, r2dot_v, flag):
    f,g = 0,0
    r2, r2dot = np.linalg.norm(r2_v), np.linalg.norm(r2dot_v)
    if (flag == "3rd" or flag == "4th"):
        u = 1/r2**3
        z = np.dot(r2_v, r2dot_v)/(r2**2)
        q = np.dot(r2dot_v, r2dot_v)/(r2**2) - u
        
        
        f = 1 - 1/2*u*tau**2 + 1/2*u*z*tau**3
        g = tau - 1/6*u*tau**3
        if (flag == "4th"):
            f += 1/24*(3*u*q-15*u*z**2+u**2)*tau**4
            g += 1/4*u*z*tau**4
        return f,g
    else:
        a = 2.2483048313568914
        a = 1 / (2/r2-r2dot**2)
        tol = flag
        n = a**(-1.5)
        delta_E, delta_E_prev = n*tau, -1
        def func(x):
            return x - (1 - r2/a)*sin(x)+(np.dot(r2_v, r2dot_v))/(n*a**2)*(1-cos(x))-n*tau
        def func_d(x):
            return 1 - (1 - r2/a)*cos(x)+(np.dot(r2_v, r2dot_v))/(n*a**2)*(sin(x))
        while (abs(delta_E - delta_E_prev) > tol):
            delta_E_prev = delta_E
            delta_E = delta_E - func(delta_E)/func_d(delta_E)
        print(delta_E)

        f,g = 1 - a/r2 * (1-cos(delta_E)), tau + 1/n*(sin(delta_E)-delta_E)
        return f,g


# test cases
tau1_1 = -0.32618569435308475
tau3_1 = 0.050840808143482484
r2_1 = [0.26640998194891174, -1.382856212643199, -0.505199925482389]
r2dot_1 = [0.8439832722802604, -0.39937767878456487, 0.14200790188593015] 

tau1_2 = -0.32618617484601165 
tau3_2 = 0.0508408854033231 
r2_2 = [0.26799552002875776, -1.3726277901924608, -0.5026729612047128]
r2dot_2 = [0.8456809141954584, -0.3838382184712308, 0.14215854191172816] 

tau1_3 = -0.3261857571141891
tau3_3 = 0.05084081855693949 
r2_3 = [0.26662393644794813, -1.381475976476564, -0.5048589337503169]
r2dot_3 = [0.8442117090940343, -0.39728396707075087, 0.14202728258915864] 

f_11, g_11 = fg(tau1_1, r2_1, r2dot_1, 1e-12)
f_31, g_31 = fg(tau3_1, r2_1, r2dot_1, 1e-12)
f_12, g_12 = fg(tau1_2, r2_2, r2dot_2, "3rd")
f_32, g_32 = fg(tau3_2, r2_2, r2dot_2, "3rd")
f_13, g_13 = fg(tau1_3, r2_3, r2dot_3, "4th")
f_33, g_33 = fg(tau3_3, r2_3, r2dot_3, "4th")

print("-----------------------------")
print("test case 1:")
print("-----------------------------")
print("f&g functions with tolerance = 1e-12 gives:")
print("f1 =", f_11)
print("f3 =", f_31)
print("g1 =", g_11)
print("g3 =", g_31)
print("-----------------------------")

print("-----------------------------")
print("test case 2:")
print("-----------------------------")
print("f&g 3rd-order series gives:")
print("f1 =", f_12)
print("f3 =", f_32)
print("g1 =", g_12)
print("g3 =", g_32)
print("-----------------------------")

print("-----------------------------")
print("test case 3:")
print("-----------------------------")
print("f&g 4rd-order series gives:")
print("f1 =", f_13)
print("f3 =", f_33)
print("g1 =", g_13)
print("g3 =", g_33)
print("-----------------------------")