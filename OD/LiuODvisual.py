import numpy as np
from math import *
from odlib import *
from vpython import *

# ast_a, ast_e, ast_i, ast_Omega, ast_omega, ast_M = a, e, radians(i), radians(Omega), radians(omega), radians(M)
ast_a, ast_e, ast_i, ast_Omega, ast_omega, ast_M = 2.739830439528865, 0.5384709236199937, 0.06905422234039398, 4.168378940765262, 1.8837245977915602, 5.939558599487618

class Body:
    name, a, e, i, Omega, omega, M, E, P = "", 0, 0, 0, 0, 0, 0, 0, 0
    pos_op, rotateO, rotatei, rotatew, pos_ec = np.zeros(3), np.zeros((3,3)), np.zeros((3,3)), np.zeros((3,3)), np.zeros(3)

    def __init__(self, name, a, e, i, Omega, omega, M, radius, color):
        self.name, self.a, self.e, self.i, self.Omega, self.omega, self.M = name, a, e, i, Omega, omega, M
        self.radius, self.color = radius, color
        self.E = nr_kepler(M, e)
        self.P = 2*pi*sqrt(a**3 / k**2)

    def calculate_pos_ec(self):
        self.pos_op = np.array([self.a*cos(self.E)-self.a*self.e, self.a*sqrt(1-self.e**2)*sin(self.E), 0])
        self.rotateO = np.array([[cos(self.Omega), -sin(self.Omega), 0], [sin(self.Omega), cos(self.Omega), 0], [0,0,1]])
        self.rotatei = np.array([[1,0,0],[0,cos(self.i), -sin(self.i)], [0,sin(self.i), cos(self.i)]])
        self.rotatew = np.array([[cos(self.omega), -sin(self.omega), 0], [sin(self.omega), cos(self.omega), 0], [0,0,1]])
        self.pos_ec = self.rotateO @ self.rotatei @ self.rotatew @ self.pos_op

    def update_pos(self, time):
        # readjust M and E based on new time
        new_M = 2*pi/self.P*time + self.M
        self.E = nr_kepler(new_M,self.e)

        #recalculate rectangular pos
        self.pos_op = np.array([self.a*cos(self.E)-self.a*self.e, self.a*sqrt(1-self.e**2)*sin(self.E), 0])

        #recalculate ecliptic coordinates
        self.pos_ec = self.rotateO @ self.rotatei @ self.rotatew @ self.pos_op
    
    def get_pos_ec(self):
        return self.pos_ec
    def get_radius(self):
        return self.radius
    def get_color(self):
        return self.color
    def get_name(self):
        return self.name

scale_1 = 1500
scale_2 = 250
def visualize():
    sun = sphere(pos=vector(0.,0.,0.), radius=0.00465040106*20, color=color.yellow)
    sun_label = label(pos = sun.pos, text = "sun")
    mercury = Body("mercury", a=0.387098, e=0.205630, i=radians(7.005), Omega=radians(48.331), omega=radians(29.124), M=ast_M, radius = 0.00001630681*scale_1, color=vector(0.888,0.736,0.540))
    venus = Body("venus", a=0.723332, e=0.0068, i=radians(3.39471), Omega=radians(76.67069), omega=radians(54.85229), M=ast_M, radius = 0.00004045454*scale_1, color=color.magenta)
    earth = Body("earth", a=1.00000261, e=0.01671123, i=0., Omega=radians(348.73936), omega=radians(114.20783), M=ast_M, radius = 0.00004263368*scale_1, color=color.blue)
    mars = Body("mars", a=1.523679, e=0.093315, i=radians(1.850), Omega=radians(49.562), omega=radians(286.537), M=ast_M, radius = 0.00002270053*scale_1, color=color.red)
    jupiter = Body("jupiter", a=5.204267, e=0.048775, i=radians(1.305), Omega=radians(100.492), omega=radians(275.066), M=ast_M, radius = 0.0004778877*scale_2, color=color.orange)
    saturn = Body("saturn", a=9.58201720, e=0.055723219, i=radians(2.485240), Omega=radians(113.642811), omega=radians(336.013862), M=ast_M, radius = 0.00040286096*scale_2, color = color.yellow)
    uranus = Body("uranus", a=19.22941195, e=0.044405586, i=radians(0.772556), Omega=radians(73.989821), omega=radians(96.541318), M=ast_M, radius = 0.00017084893*scale_2, color = color.cyan)
    neptune = Body("neptune", a=30.10366151, e=0.011214269, i=radians(1.767975), Omega=radians(131.794310), omega=radians(265.646853), M=ast_M, radius = 0.00016553475*scale_2, color = color.blue)
    asteroid = Body("2003 UD8", a=ast_a, e=ast_e, i=ast_i, Omega=ast_Omega, omega=ast_omega, M=ast_M, radius = 0.00001*scale_1, color = color.white)
    bodies = [[mercury], [venus], [earth], [mars], [jupiter], [saturn], [uranus], [neptune], [asteroid]]

    # in the array "bodies", each element is itself an array:
    # 0th element = body object containing all orbital element information
    # 1st element = sphere object in vpython 
    # 2nd element = label object in vpython
    for i in range(len(bodies)):
        bodies[i][0].calculate_pos_ec()
        bodies[i].append(sphere(pos=vector(bodies[i][0].get_pos_ec()[0], bodies[i][0].get_pos_ec()[1], bodies[i][0].get_pos_ec()[2]), radius=bodies[i][0].get_radius(), color=bodies[i][0].get_color()))
        bodies[i][1].trail = curve(color = bodies[i][0].get_color())
        bodies[i].append(label(pos=bodies[i][1].pos, text=bodies[i][0].get_name()))

    time = 0
    print("------------------")
    print("  ^~^  ,")
    print(" ('Y') )")
    print(" /   \/ ")
    print("(\|||/) ")
    print("------------------")
    print("~thanks for running!")

    while (True):
            rate(100)
            time = time + 1

            for i in range(len(bodies)):
                bodies[i][0].update_pos(time)
                bodies[i][1].pos = vector(bodies[i][0].get_pos_ec()[0], bodies[i][0].get_pos_ec()[1], bodies[i][0].get_pos_ec()[2])
                bodies[i][2].pos = bodies[i][1].pos
                bodies[i][1].trail.append(pos=bodies[i][1].pos)  