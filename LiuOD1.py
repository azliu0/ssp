import numpy as np

from odlib import k, retrieve_data, ang_momentum

#extracting final answer
pos,velocity = retrieve_data("LiuInput.txt", "2018-Jul-14", "00:00:00.0000")

print([round(element, 6) for element in ang_momentum(pos, velocity)])