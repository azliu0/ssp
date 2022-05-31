import numpy as np


circle, points, x , y, count = 0, 0 , 0 , 0, 0
pi_list = []

while(True):
    while (points < 1e6):
        x,y = 2*np.random.random(), 2*np.random.random()
        if ((x - 1)**2 + (y-1)**2 <= 1):
            circle += 1
        points += 1
    pi_list.append(4*circle/points)
    count += 1
    print("mean:", np.sum(pi_list)/count)
    points, circle = 0,0

