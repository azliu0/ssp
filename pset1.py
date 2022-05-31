# SSP 2021 Python pset 1
# 06/23/21
# Andrew Liu

from math import *

# problem 1
def min3(a,b,c): 
    values = [a,b,c]
    values.sort()
    return values[0]
def max3(a,b,c):
    values = [a,b,c]
    values.sort()
    return values[2]
def spread3(a,b,c):
    return max3(a,b,c)-min3(a,b,c)


# problem 2
def convert_angle(degrees, minutes, seconds, rads, normalize):
    # handle a negative angle
    minutes = copysign(minutes,degrees)
    seconds = copysign(seconds, degrees)

    # perform angle conversion
    angle_degrees = degrees + minutes/60 + seconds/3600

    if normalize:
        while(angle_degrees>= 360):
            angle_degrees = angle_degrees - 360
        while(angle_degrees < 0):
            angle_degrees = angle_degrees+360



    # return result
    if rads:
        return radians(angle_degrees)
    return angle_degrees


# problem 3
# RA and dec should be in decimal degrees
# L and lat are longitude and latitude of observer in decimal degrees
# Y = 4-digit year; 1901 ≤ Y ≤ 2099
# M = month; 1 ≤ M ≤ 12
# D = day; 1 ≤ D ≤ 31
# UT should be a UT decimal hour
def RADec_to_AltAz(RA, dec, lon, lat, Y, M, D, UT):
    
    #calculating julian date @ UT=0:
    J = 367*Y - (7*(Y+(M+9)//12))//4 + (275*M)//9+D+1721013.5
    
    #calculating LST:
    J = (J-2451545)/36525
    GMST = 100.46061837 + 36000.77053608*J+3.87933*(10**(-4))*(J**2)-(J**3)/(3.871*10**7)
    LST = GMST + 360.985647366*(UT/24)+lon

    #calculating hour angle (h): 
    h = LST-RA #hour angle should be a value in the range (-180,180]
    while (h > 180):
        h = h-360
    while (h <= -180):
        h = h+360

    #calculating altitude (al):
    #no need to check for ranges because abs(al) <= 90 always
    al = asin(cos(radians(h))*cos(radians(lat))*cos(radians(dec)) + sin(radians(lat))*sin(radians(dec))) 

    #calculating azimuth (az):
    az = asin(-cos(radians(dec))*sin(radians(h))/cos(al))

    al = degrees(al)
    az = degrees(az)

    # we need to adjust azimuth because range of possibilities for az is (0,360), which superseceds the domain of the inverse sine function
    #
    # calculating the length of arc from "N" to the star, by using the cosine formula on the spherical triangle with vertices 
    # N, NCP, and our star.
    #
    # if the length is >90, then need to return 180-az, because the star is in the southern sky
    # otherwise, we take az as it is, because it is in the northern sky (ie the primary branch of the domain of the inverse sin function).

    arc = degrees(acos(cos(radians(lat))*sin(radians(dec)) - cos(radians(h))*sin(radians(lat))*cos(radians(dec))))

    if arc > 90:
            return al,180-az
    else:
        if (az<0):
            return al, az+360
        return al,az



# lats and lons are in decimal degrees, negative indicates west/south
def great_circle_dist(lat1, lon1, lat2, lon2):
    R=6371000
    a = sin(radians((lat2-lat1)/2))**2+cos(radians(lat1))*cos(radians(lat2))*sin(radians((lon2-lon1)/2))**2
    c=2*atan2(sqrt(a), sqrt(1-a))
    d=R*c
    return d  