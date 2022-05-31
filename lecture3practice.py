import numpy as np
import math

fruits = np.array([["Apple", "Banana", "Blueberry", "Cherry"],
                   ["Coconut", "Grapefruit", "Kumquat", "Mango"],
                   ["Nectarine", "Orange", "Tangerine", "Pomegranate"],
                   ["Lemon", "Raspberry", "Strawberry", "Tomato"]])

print(fruits[3,3])
print(fruits[1:3, 1:3])
# print(fruits[0:3:2])
# print(fruits[2:0:-1, 2:0:-1])
print(fruits[:,::-1])
# copyfruits = fruits.copy()
print(fruits[1].shape[0])

copy = fruits[1:3,1:3]
print(copy)
print(fruits)
# fruits = np.array(["Apple", "Banana", "Blueberry", "Cherry"])
# print(fruits.shape)
# print(copyfruits)
