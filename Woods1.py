# -*- coding: utf-8 -*-
"""
Created on Sun Jun 11 20:32:12 2017

@author: Daniel
"""

import numpy                 #loading our favorite library
from matplotlib import pyplot    #and the useful plotting library
from matplotlib import animation

height = 2
length = 1
########################################################
def init():
    line.set_data([], [])
    return line,

fig = pyplot.figure()
ax = pyplot.axes(xlim=(0, length), ylim=(0, height))
line, = ax.plot([], [], lw=2)
#########################################################
    

nx = 101
dx = length / float(nx - 1)
nt = 20    #the number of timesteps we want to calculate
nu = 0.3   #the value of viscosity
Pe = 6993006
lamb = 1.2
#dt = sigma * dx**2 / nu #dt is defined using sigma ... more later!
dt = 0.005

T = numpy.zeros(nx)      #a numpy array with nx elements all equal to 1.
#u[int(0.5 / dx):int(1.0 / dx + 1)] = 2  #setting u = 2 between 0.5 and 1 as per our I.C.s
T[0] = 1

Tn = numpy.ones(nx) #our placeholder array, un, to advance the solution in time



#for n in range(nt):  #iterate through time
def animate(j):
    #T[0:int(0.1*j)]=1
    
    Tn = T.copy() ##copy the existing values of u into un
    for i in range(1, nx - 1):
    #for i in range(1, nx - 1):
        T[i] = Tn[i] - lamb * (dt /dx) * (Tn[i] - Tn[i-1]) +  (dt / dx**2) * (1/Pe) * (Tn[i+1] - 2 * Tn[i] + Tn[i-1])

        line.set_data(numpy.linspace(0, length, nx), T)
        
    return line,
    
    #pyplot.plot(numpy.linspace(0, 2, nx), u);
    
# call the animator.  blit=True means only re-draw the parts that have changed.
anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=20000, interval=20, blit=True)