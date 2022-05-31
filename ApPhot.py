# Programming Assignment #5 / SSP
# Andrew Liu
# Tuesday, July 6th, 2021

from math import sqrt, ceil, floor, log
import numpy as np
from astropy.io import fits

gain = 0.8
dark = 10
read = 11

# find sum and number of non-zero elements in array 
def sumSize(array):
    sum, size = 0,0
    for row in array:
        for element in row:
            if (element):
                sum += element
                size += 1
    return sum,size

# crop array from center and radius according to center-pixel inclusion test 
def centerCrop(inputArray, center_x, center_y, radius):
    arrayCopy = inputArray.copy()
    arrayCopy = arrayCopy[floor(center_y)-radius:ceil(center_y)+radius+1, floor(center_x)-radius:ceil(center_x)+radius+1]
    for y in range(len(arrayCopy)):
        for x in range(len(arrayCopy[0])):
            if ((x-(center_x - floor(center_x) + radius))**2 + (y-(center_y - floor(center_y) + radius))**2 > radius**2):
                arrayCopy[y,x]=0
    return arrayCopy

# centroid calculation function
def centroid(array):
    weighted_sum, total_sum = 0,0

    #calculating the centroid location, centroid
    for row in array:
        for index, element in enumerate(row):
            weighted_sum = weighted_sum + (element)*index
            total_sum = total_sum + (element)
    centroid = weighted_sum / total_sum

    return centroid

# finding the centroid of target object, input array is the raw fits file array
def findCentroid(inputArray, target_x, target_y, radius):

    #crop array into centered array of interest
    arrayCopy = centerCrop(inputArray, target_x, target_y, radius)

    centroid_x = centroid(arrayCopy)
    centroid_y = centroid(arrayCopy.T)

    return centroid_x+target_x-radius, centroid_y+target_y-radius

# main method
def apPhot(fits_file, target_x, target_y, ap, i_an, o_an):
    arrayPix = np.array(fits.getdata(fits_file))

    #centroid values
    centroid_x,centroid_y = findCentroid(arrayPix, target_x, target_y, ap+2)

    #calculating ADU and n values for aperture 
    aperture = centerCrop(arrayPix, centroid_x, centroid_y, ap)
    ADU_ap, n_ap = sumSize(aperture)

    #calculating ADU and n values for annulus. need extra steps because size of inner annulus != size of outer annulus
    inner_annulus = centerCrop(arrayPix, centroid_x, centroid_y, i_an)
    outer_annulus = centerCrop(arrayPix, centroid_x, centroid_y, o_an)
    pad = o_an - i_an
    temp = np.zeros(outer_annulus.shape)
    temp[pad:inner_annulus.shape[0]+pad, pad:inner_annulus.shape[1]+pad] = inner_annulus
    annulus = outer_annulus-temp
    ADU_an, n_an = sumSize(annulus)

    # calculating final values
    S_ADU = ADU_ap - ADU_an*(n_ap / n_an)
    m_inst = -2.5*log(S_ADU, 10)

    print(ADU_ap, ADU_an, n_ap, n_an)
    print(S_ADU, gain, dark, read)
    SNR = sqrt(S_ADU*gain/(1+n_ap*(1+n_ap/n_an)*(dark+ADU_an/n_an*gain + read**2+gain**2/12)/(S_ADU*gain)))
    delta_m = 1.0875/SNR

    print("-------------------------------------------")
    print("centroid = ", centroid_x, ",", centroid_y)
    print("Signal = ", S_ADU, "; SNR = ", SNR)
    print("inst mag =", m_inst, "+/-", delta_m)
    print("-------------------------------------------")

apPhot("aptest.FIT", 93, 285, 3,5,10)
apPhot("aptest.FIT", 490, 293, 5,8,13)
