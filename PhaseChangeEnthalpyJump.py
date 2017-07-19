# -*- coding: utf-8 -*-
"""
Created on Sun Jun 11 20:32:12 2017

@author: Daniel
"""

import numpy                 #loading our favorite library
from matplotlib import pyplot    #and the useful plotting library
from matplotlib import animation
import time


t0 = time.time()
height = 1
length = 1
########################################################
def init():
    line.set_data([], [])
    return line,

#fig = pyplot.figure()
#ax = pyplot.axes(xlim=(0, length), ylim=(0, height))
#line, = ax.plot([], [], lw=2)
#########################################################
    
for w in range(1,2):
    nx = 501                   # Number of mesh points
    dx = length / float(nx - 1) # Distance between mesh points
    nt = 2400                     # The number of timesteps we want to calculate
    dt = 0.0001                   # Size of the time steps
    
    Stg = 1000.0                  # Stanton Number of the gas phase
    Stl = 200.0                   # Stanton Numer of the liquid phase
    Sts = 4000.0                   # Stanton Numer of the solid phase
    Stjump = 80.0                # Stanton Numer of the jump conditional phase
    rho_lv = 141.0                # Ratio of the liquid phase to the vapor phase
    Tsat = 0.99
    delta_theta = 1.0 - Tsat    # Change in theta in the jump conditional phase e.g. (1-0.99)
    rho = 141.0                   # Non-dimensionalized density ratio from liquid to vapor
    vel = 1.0                     # Non-dim. velocity
    
    Tl = numpy.zeros(nx)     # Initializing temperature of the liquid phase
    Ts = numpy.zeros(nx)     # Initializing temperature of the solid phase
    f = numpy.ones(nx)       # Initializing vapor quality
    rho = numpy.ones(nx)     # Initializing density
    dTdt = numpy.zeros(nx)   # Initializing derivative of temperature of with time
    v = numpy.zeros(nx)       # Initializing velocity
    dvdx = numpy.zeros(nx)   # Initializing derivative of velocity with distance
    
    v[0] = 1
    Tl[0] = 1
    Ts[0] = 1
    #u[int(0.5 / dx):int(1.0 / dx + 1)] = 2  #setting u = 2 between 0.5 and 1 as per our I.C.s

    
    Tln = numpy.ones(nx) #our placeholder array, Tln, to advance the solution in time
    Tsn = numpy.ones(nx) #our placeholder array, Tsn, to advance the solution in time

    for n in range(nt):  #iterate through time
    #def animate(j):
        Tln = Tl.copy()
        Tsn = Ts.copy() ##copy the existing values of T into Tn
        
        for i in range(1, nx - 1):
            Tl[i] = Tln[i] - Stg * dt * (Tln[i]-Tsn[i]) - (dt/float(dx))*(Tln[i]-Tln[i-1])
            Ts[i] = Tsn[i] - Sts * dt * (Tln[i]-Tsn[i])
            
            if Tl[i]>1.0:
                rho[i] = 1.0
                v[i] = v[i-1]
            
            if Tl[i]<1.0 and Tl[i]>Tsat:
                rho[i] = 1 + (1-(Tl[i]-Tsat)/delta_theta)*rho_lv
                dTdt[i] = -Stjump*(Tl[i]-Ts[i])-(Tl[i]-Tl[i-1])/dx
                v[i] = v[i-1] + (((rho_lv/delta_theta)*(dTdt[i]+v[i-1]*(Tl[i]-Tl[i-1])/dx))/rho[i])*dx
                dvdx[i] = (((rho_lv/delta_theta)*(dTdt[i]+v[i]*(Tl[i]-Tl[i-1])/dx))/rho[i])
                f[i] = (rho_lv/delta_theta)*(Tl[i]-Tsat)*(dTdt[i]+v[i]*((Tl[i]-Tl[i-1])/dx))
                Tl[i] = Tl[i] + dt*((Stjump*(-Tl[i]+Ts[i]))/rho[i]+f[i]/rho[i]-(v[i]*(Tl[i]-Tl[i-1])/dx)-(Tl[i]-Tsat)*dvdx[i])
                Ts[i] = Ts[i-1] + Sts * dt * (Tl[i]-Ts[i])
            if Ts[i] < 0.99:
                rho[i] = 141.0
                v[i] = v[i-1]*rho[i-1]/rho[i]
                Tl[i] = Tl[i] - Stl * dt * (Tl[i]-Ts[i]) - (dt/dx)*(Tl[i]-Tl[i-1])
                Ts[i] = Ts[i] + Sts * dt * (Tl[i]-Ts[i])
            v[nx-1] = v[nx-2]
            

                
    #        line.set_data(numpy.linspace(0, length, nx), theta)
    #        
    #    return line,
        #pyplot.figure(1)
        #pyplot.semilogx(numpy.linspace(0, 2, nx), theta)
#    pyplot.figure(2)
#    pyplot.semilogx(numpy.linspace(0, length, nx), Tl)
#    pyplot.semilogx(numpy.linspace(0.001, length, nx), P)
    pyplot.plot(numpy.linspace(0, 1, nx), Tl);
    pyplot.plot(numpy.linspace(0, 1, nx), v);
        
    # call the animator.  blit=True means only re-draw the parts that have changed.
    #anim = animation.FuncAnimation(fig, animate, init_func=init,
    #                               frames=20000, interval=20, blit=True)
t1 = time.time()
total = t1-t0
print total