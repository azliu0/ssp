# SSP 2021 Python pset 2
# Andrew Liu
# 6-26-2021

from math import sin, cos, asin, acos, pi, radians, degrees

# a function to determine the quadrant of an angle based on its sine and cosine (in radians)
# returns the angle in the correct quadrant (in radians)
def findQuadrant(sine, cosine):
    if cosine > 0 and sine > 0: #1
        return asin(sine)

    if cosine < 0 and sine > 0: #2
        return acos(cosine)

    if cosine < 0 and sine < 0: #3
        return pi - asin(sine)

    if cosine > 0 and sine < 0: #4
        return 2*pi + asin(sine)

# a function that given the values (in decimal degrees) of two sides and the included angle of a spheical triangle
# returns the values of the remaining side and two angles (in decimal degrees)
def SAS(a, B, c):
    a = radians(a)
    B = radians(B)
    c = radians(c)

    b = acos(cos(B)*sin(c)*sin(a) + cos(c)*cos(a))
    A = acos((cos(a)-cos(c)*cos(b)) / (sin(c) * sin(b)))
    C = acos((cos(c)-cos(a)*cos(b)) / (sin(a) * sin(b)))
    return degrees(b), degrees(A), degrees(C) 

# DO NOT REMOVE OR MODIFY THIS CODE
s3, a1, a2 = SAS(106, 114, 42)
if abs(s3 - 117.804) > 1e-3 or abs(a1 - 83.11) > 1e-3 or abs(a2 - 43.715) > 1e-3:
    print("SAS function INCORRECT, expected (117.804, 83.11, 43.715), but got", (s3, a1, a2))
else:
    print("SAS function CORRECT")
