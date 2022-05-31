# Programming Assignment #4 / SSP
# Andrew Liu
# Thursday, July 1st, 2021
# NOTE: this code passed both test cases exactly as intended

import numpy as np
from math import sqrt

# ------------------------------------------------------------------------------
# helper functions for the main method:
# ------------------------------------------------------------------------------

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

#functions for reversing the above two functions (convert degrees to sg)
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

# this function takes two arrays and returns the sum of the products of each of their corresponding elements
def combine(arr1, arr2):
    sum = 0
    for i in range(len(arr1)):
        sum += arr1[i]*arr2[i]
    return sum

# ------------------------------------------------------------------------------
# main method:
# ------------------------------------------------------------------------------

def LSPR(N, file, x_test, y_test):
    #read in file
    test_case = open(file, "r")

    #initialize arrays for x pixel, y pixel, right ascensions, and declinations
    x_pix, y_pix, ra, dec = [],[],[],[]

    # fully read in the data, and sort into respective arrays
    i=0
    for line in test_case:
        for word in line.split():
            if (i%4 == 0):
                x_pix.append(float(word))
            elif (i%4 == 1):
                y_pix.append(float(word))
            elif (i%4 == 2):
                ra.append(ra_to_degrees(word))
            elif (i%4 == 3):
                dec.append(dec_to_degrees(word))
            i += 1
    x_pix,y_pix,ra,dec=np.array(x_pix),np.array(y_pix),np.array(ra),np.array(dec)


    #matrix used to calculate RA and dec: 
    ones = np.ones(N)
    eqn_matrix = np.array([
        [N, combine(x_pix, ones), combine(y_pix, ones)],
        [combine(x_pix, ones), combine(x_pix, x_pix), combine(x_pix, y_pix)],
        [combine(y_pix, ones), combine(x_pix, y_pix), combine(y_pix, y_pix)]
    ])

    #calculating RA plate constants:
    ra_matrix = np.array([combine(ra, ones), combine(ra, x_pix), combine(ra, y_pix)]).reshape(3,1)
    ra_plate = np.linalg.inv(eqn_matrix) @ ra_matrix 

    #calculating dec plate constants:
    dec_matrix = np.array([combine(dec, ones), combine(dec, x_pix), combine(dec, y_pix)]).reshape(3,1)
    dec_plate = np.linalg.inv(eqn_matrix) @ dec_matrix 

    #calculating standard deviation for ra and dec
    delta_ra = 0
    delta_dec = 0
    for i in range(N):
        delta_ra += (ra[i]-ra_plate[0,0]-ra_plate[1,0]*x_pix[i]-ra_plate[2,0]*y_pix[i])**2
        delta_dec += (dec[i]-dec_plate[0,0]-dec_plate[1,0]*x_pix[i]-dec_plate[2,0]*y_pix[i])**2
    delta_ra = sqrt(1/(N-3) * delta_ra) * 3600
    delta_dec = sqrt(1/(N-3) * delta_dec) * 3600

    #calculating and formatting final ra,dec for given x,y
    ra_final = ra_plate[0,0]+ra_plate[1,0]*x_test+ra_plate[2,0]*y_test
    dec_final = dec_plate[0,0] + dec_plate[1,0]*x_test + dec_plate[2,0]*y_test

    #formatting print statements:
#     print('***************')
#     print('plate constants')
#     print('***************')
#     print('b1:', ra_plate[0,0], 'deg')
#     print('b2:', dec_plate[0,0], 'deg')
#     print('a11:', ra_plate[1,0], 'deg/pix')
#     print('a12:', ra_plate[2,0], 'deg/pix')
#     print('a21:', dec_plate[1,0], 'deg/pix')
#     print('a22:', dec_plate[2,0], 'deg/pix')
#     print('***********')
#     print('uncertainty')
#     print('***********')
#     print('RA:', "{:.2f}".format(delta_ra), 'arcsec')
#     print('dec:', "{:.2f}".format(delta_dec), 'arcsec')
#     print('***********************************')
#     print('astrometry for')
#     print('(x,y)=', "("+str(x_test)+','+str(y_test)+")")
#     print('***********************************')
#     print('RA = ', ra_to_sg(ra_final))
#     print('dec = ', dec_to_sg(dec_final))
    


# print('------------------------------------------------------------------------------')
# print('first test case')
# print('------------------------------------------------------------------------------')
# LSPR(12, "LSPRtestinput1.txt", 484.35, 382.62)
# print('------------------------------------------------------------------------------')
# print('second test case')
# print('------------------------------------------------------------------------------')
# LSPR(12, "LSPRtestinput2.txt", 1403.6, 1585.9)


