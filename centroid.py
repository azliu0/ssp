# Programming Assignment #3 / SSP
# Andrew Liu
# Tuesday, June 29th, 2021

from math import sqrt
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from vpython import *

#this method calculates the x_centroid and x_delta for any arbitrary array
def centroidAndDelta(array):
    weighted_sum, delta, total_sum = 0,0,0

    #calculating the centroid location, centroid
    for row in array:
        for index, element in enumerate(row):
            weighted_sum = weighted_sum + element*index
            total_sum = total_sum + element
    centroid = weighted_sum / total_sum

    #calculating the standard deviation, delta
    for row in array:
        for index, element in enumerate(row):
            delta = delta + element*((index - centroid)/total_sum)**2
    delta = sqrt(delta)

    return centroid,delta

#this method calculates the x,y centroids and standard deviations for a given pix file, target, and radius
def findCentroid(fits_file, target_x, target_y, radius):
    arrayPix = np.array(fits.getdata(fits_file))
    
    #crop array into centered array of interest
    arrayPix = arrayPix[target_y-radius:target_y+radius+1, target_x-radius:target_x+radius+1]

    centroid_x, delta_x = centroidAndDelta(arrayPix)
    centroid_y, delta_y = centroidAndDelta(arrayPix.T)

    return centroid_x+target_x-radius, centroid_y+target_y-radius, delta_x, delta_y 

#-------------------------------------------------------------------------------------------------------------------------------------
#test case code:
#-------------------------------------------------------------------------------------------------------------------------------------

centroid_x, centroid_y, uncert_x, uncert_y = findCentroid("sampleimage.fits", 351, 154, 1)
if abs(centroid_x - 350.9958) < 1e-3 and abs(centroid_y - 153.9955) < 1e-3:
    print("centroid calculation CORRECT")
else:
    print(f"centroid calculation INCORRECT, expected (350.9958, 153.9955), got ({centroid_x}, {centroid_y})")
if abs(uncert_x - 0.005254018) < 1e-6 and abs(uncert_y - 0.005249733) < 1e-6:
    print("uncertainty calculation CORRECT")
else:
    print(f"uncertainty calculation INCORRECT, expected (0.005254018, 0.005249733), got ({uncert_x}, {uncert_y})")

#-------------------------------------------------------------------------------------------------------------------------------------
#visualization: 
#-------------------------------------------------------------------------------------------------------------------------------------

a,b,c,d = findCentroid("aptest.FIT", 93, 285, 5)
print(a,b)

# target_x = 351
# target_y = 154
# radius=1

# centroid_x = centroid_x + radius - target_x
# centroid_y = centroid_y + radius - target_y

# arrayPix = np.array(fits.getdata("sampleimage.fits"))
# arrayPix = arrayPix[target_y-radius:target_y+radius+1, target_x-radius:target_x+radius+1]

# #calculating mean value of every element
# total_sum = 0
# for row in arrayPix:
#     for element in row:
#         total_sum = total_sum + element
# mean_radius = total_sum / arrayPix.size

# #printing array of spheres
# sphere(pos = vector(centroid_x, centroid_y,0.015), radius=mean_radius*6e-5, color = color.red)

# for i in range(arrayPix.shape[0]):
#     for j in range(arrayPix.shape[1]):
#         sphere(pos = vector(i,j,0), radius=arrayPix[i,j]*6e-5, color = color.yellow)



