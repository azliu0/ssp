'''
tuple: immutable sequence
list: mutable sequence
'''

import vpython

vpython.scene.autoscale = False

ball = vpython.sphere(radius=0.5)
ball.velocity = vpython.vector(25,5,0)

wall_r = vpython.box(pos = vpython.vector(6,0,0), size = vpython.vector(0.2,12,12), color = vpython.color.red)
wall_l = vpython.box(pos = vpython.vector(-6,0,0), size = vpython.vector(0.2,12,12), color = vpython.color.green)

vel_arr = vpython.arrow(pos = ball.pos, axis = ball.velocity, color = vpython.color.yellow)

delta_t = 0.005
while True:
    vpython.rate(50)
    ball.pos = ball.pos + ball.velocity * delta_t

    if ball.pos.x + ball.radius> wall_r.pos.x:
        ball.velocity.x = -ball.velocity.x
    if ball.pos.x - ball.radius < wall_l.pos.x:
        ball.velocity.x = -ball.velocity.x

    vel_arr.pos = ball.pos
    vel_arr.axis = ball.velocity * 0.2