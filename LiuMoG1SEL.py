# Andrew Liu MoG1: SEL
# SSP July 9th, 2021

import numpy as np

# NOTE: code passes both test cases perfectly

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

taus_1 = [-0.15481889055,0.15481889055, 0.3096377811] 
sun2_1 = [-0.2398478458274071, 0.9065739917845802, 0.3929623749770952] 
rhohat2_1 = [-0.8518563498182248, -0.2484702599212149, 0.4610892421311239]
Ds_1 = [-0.0010461861084885213, -0.17297581974209159, -0.17201260125558127, -0.16712421570714076]

taus_2 = [-0.1720209895, 0.1720209895, 0.344041979] 
sun2_2 = [-0.2683394448727136, 0.8997620807182745, 0.3900022331276332] 
rhohat2_2 = [0.052719013914983195, -0.9121555187306237, 0.40643943610469035] 
Ds_2 = [0.0011797570297704812, 0.052586743761143424, 0.05848153743706686, 0.06274019190783499] 

roots_1, rhos_1 = SEL(taus_1, sun2_1, rhohat2_1, Ds_1)
roots_2, rhos_2 = SEL(taus_2, sun2_2, rhohat2_2, Ds_2)

print("-------------------------------------------------")
print("test case 1:")
print("-------------------------------------------------")
print("roots:", roots_1)
print("rhos:", rhos_1)

print("-------------------------------------------------")
print("test case 2:")
print("-------------------------------------------------")
print("roots:", roots_2)
print("rhos:", rhos_2)


