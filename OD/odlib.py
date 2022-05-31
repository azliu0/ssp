import numpy as np
from math import pi, sqrt, atan, asin, acos, sin, cos, degrees, radians
from astropy.io import fits

#gaussian constant = regular days to gauss days conversion
k = 0.0172020989484
cAU = 173.144643267 # in units of AU / day
eps = radians(23.4373)

#retrieve data for position and velocity vectors for given jpl file, date, and time
def retrieve_data(file, date, time):

    jpl_data = open(file, "r")
    data=[]

    # read in the data, up until the desired date
    for line in jpl_data:
        if (str(date) + " " + str(time)) in line:
            break

    # we want the next two lines, which will give us our position and velocity vectors
    data.append(jpl_data.readline()[:-1].split(" "))
    data.append(jpl_data.readline()[:-1].split(" "))

    pos = np.array([float(data[0][3]), float(data[0][6]), float(data[0][9])])
    velocity = np.array([float(data[1][2]), float(data[1][4]), float(data[1][6])])*1/k
    
    return pos,velocity

# this functions returns an angle given its sin and cosine value
def quadrant_check(sin_v, cos_v):
    tan_v = sin_v / cos_v
    theta = atan(tan_v)
    if (cos_v < 0):
        while(theta < pi/2):
            theta += pi
        return theta
    else:
        while(theta < 0):
            theta += 2*pi
        return theta

def juliandate(Y, M, D, UT):
    return 367*Y - (7*(Y+(M+9)//12))//4 + (275*M)//9+D+1721013.5 + UT/24

#functions for reading in values (convert sg to degrees)
def ra_to_degrees(word):
    sg = word.split(":")
    return float(sg[0])*15+float(sg[1])*15/60+float(sg[2])*15/3600

def dec_to_degrees(word):
    sign = True
    if word[0] == "-":
        sign = False
    word = word[1:]
    sg = word.split(":")
    ret_value = float(sg[0])+float(sg[1])/60+float(sg[2])/3600
    if sign:
        return ret_value
    else:
        return (-1)*ret_value  

def ra_to_sg(value):
    value = value / 15
    hours = value // 1
    value = (value - hours)*60
    minutes = value // 1
    value = (value - minutes)*60
    seconds = value
    return "{:.0f}".format(hours)+':'+"{:.0f}".format(minutes)+':'+"{:.2f}".format(seconds)

def dec_to_sg(value):
    sign = True
    if(value < 0):
        sign = False
        value = (-1)*value
    degrees = value // 1
    value = (value - degrees)*60
    arcminutes = value // 1
    value = (value - arcminutes)*60
    arcseconds = value
    ret_string = "{:.0f}".format(degrees)+':'+"{:.0f}".format(arcminutes)+':'+"{:.1f}".format(arcseconds)
    if sign:
        return '+'+ret_string
    else:
        return '-'+ret_string

# function for converting equatorial spherical coordinates into equatorical cartesian coordinates
def sph_to_cart(ra, dec):
    z = sin(dec)
    x = cos(ra)*cos(dec)
    y = sin(ra)*cos(dec)
    return np.array([x,y,z])

#h is defined as r cross r dot
def ang_momentum(pos, velocity):
    return np.cross(pos, velocity)

# calculates and returns the six orbital elements: 
# semimajor axis, eccentricity, inclination, long of ascending node, argument of perihelion, mean anomaly
def orbital_elements(pos,velocity):
    am = ang_momentum(pos,velocity)
    h = sqrt(np.sum(am*am))
    r = sqrt(np.sum(pos*pos))
    v = sqrt(np.sum(velocity*velocity))

    #semimajor_axis
    a = 1 / (2/r-v**2)

    #eccentricity
    e = sqrt(1 - h**2/a)

    #inclination
    i = atan(sqrt(am[0]**2+am[1]**2) / am[2])

    #longitude of ascending node
    sin_long_an, cos_long_an = am[0]/(h*sin(i)), -am[1]/(h*sin(i))
    Omega = quadrant_check(sin_long_an, cos_long_an)

    #argument of perihelion and true anomaly
    sin_U,cos_U = pos[2]/(r*sin(i)), (pos[0]*cos(Omega) + pos[1]*sin(Omega))/r
    U = quadrant_check(sin_U, cos_U)

    #calculating nu_0
    sin_ta, cos_ta = a*(1-e**2)/(e*h) * (np.sum(pos * velocity))/r, 1/e * (a*(1-e**2)/r - 1)
    nu_0 = quadrant_check(sin_ta, cos_ta)

    omega = U - nu_0
    while (omega < 0):
        omega += 2*pi
    while (omega > 2*pi):
        omega -= 2*pi

    #eccentric anomaly
    E = acos(1/e * (1 - r/a))
    if nu_0 > pi:
        E = 2*pi - E
    
    #mean anomaly
    M = E - e*sin(E)
    while (M < 0):
        M += 2*pi
    while (M > 2*pi):
        M -= 2*pi
    
    return a, e, degrees(i), degrees(Omega), degrees(omega), degrees(M)

# same exact function as above, except returns eccentric anomaly as well
def orbital_elements_with_E(pos,velocity):
    am = ang_momentum(pos,velocity)
    h = sqrt(np.sum(am*am))
    r = sqrt(np.sum(pos*pos))
    v = sqrt(np.sum(velocity*velocity))

    #semimajor_axis
    a = 1 / (2/r-v**2)

    #eccentricity
    e = sqrt(1 - h**2/a)

    #inclination
    i = atan(sqrt(am[0]**2+am[1]**2) / am[2])

    #longitude of ascending node
    sin_long_an, cos_long_an = am[0]/(h*sin(i)), -am[1]/(h*sin(i))
    Omega = quadrant_check(sin_long_an, cos_long_an)

    #argument of perihelion and true anomaly
    sin_U,cos_U = pos[2]/(r*sin(i)), (pos[0]*cos(Omega) + pos[1]*sin(Omega))/r
    U = quadrant_check(sin_U, cos_U)



    #calculating nu_0
    sin_ta, cos_ta = a*(1-e**2)/(e*h) * (np.sum(pos * velocity))/r, 1/e * (a*(1-e**2)/r - 1)
    nu_0 = quadrant_check(sin_ta, cos_ta)

    omega = U - nu_0
    while (omega < 0):
        omega += 2*pi
    while (omega > 2*pi):
        omega -= 2*pi

    #eccentric anomaly
    E = acos(1/e * (1 - r/a))
    if nu_0 > pi:
        E = 2*pi - E
    
    #mean anomaly
    M = E - e*sin(E)
    while (M < 0):
        M += 2*pi
    while (M > 2*pi):
        M -= 2*pi
    
    return a, e, degrees(i), degrees(Omega), degrees(omega), degrees(M), degrees(E)


def nr_kepler(M, e):
    tol = 0.5e-12
    E_guess, E_guess_prev = M,-1
    while (abs(E_guess - E_guess_prev) > tol):
        E_guess_prev = E_guess
        E_guess = E_guess - (M - (E_guess - e*sin(E_guess)))/(e*cos(E_guess) - 1)
    return E_guess

# returns whether or not x is an integer, used for reading in user input data
def is_int(x):
        try:
            x = int(x)
        except ValueError:
            return False
        return True

'''
calculating RMS of differences between indexed RA/dec stars and actual field values. use resulting value to generate normal distribution
of error values and simulate differences in orbital elements.
'''

def RMS_uncertainty(file):
    table = fits.open(file)[1].data
    index_ra, index_dec, field_ra, field_dec = np.array(table.index_ra), np.array(table.index_dec), np.array(table.field_ra), np.array(table.field_dec)

    ra_diff = index_ra - field_ra
    dec_diff = index_dec - field_dec

    ra_sigma = sqrt(np.mean(ra_diff**2))
    dec_sigma = sqrt(np.mean(dec_diff**2))

    return radians(ra_sigma), radians(dec_sigma)

def error(actual, calculated):
    return 100*abs(actual - calculated)/actual




    
