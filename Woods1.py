# -*- coding: utf-8 -*-
"""
Created on Sun Jun 11 20:32:12 2017

@author: Daniel
"""

import numpy                 #loading our favorite library
from matplotlib import pyplot    #and the useful plotting library
from matplotlib import animation

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
    nx = 10001
    dx = length / float(nx - 1)
    nt = 10000   #the number of timesteps we want to calculate
    Pe = 10.0
    Pe2 = 0.1 
    lamb = 1.37
    dt = 0.000001
    Length = 100.0
    kappa = 0.143E-6
    Q = (0.143*0.006)*10**-6
    phi = 0.01
    rho_l = 1000.0 #Density of liquid
    L = 2260.0     #Latent Heat of vaporization
    rho_p = 800.0
    C_p = 3000.0
    f = 0.4
    k = 10E-13
   

    #dx = numpy.logspace(-3,2,nx)
    #dt = numpy.logspace(-7,1,nx)
    
    T = numpy.ones(nx)*2      #a numpy array with nx elements all equal to 1.
    theta = numpy.ones(nx)*2
    theta[0] = 0
    #u[int(0.5 / dx):int(1.0 / dx + 1)] = 2  #setting u = 2 between 0.5 and 1 as per our I.C.s
    T[0] = 0.01
    
    Tn = numpy.ones(nx) #our placeholder array, Tn, to advance the solution in time
    
    
    P = numpy.ones(nx)
    Pn = numpy.ones(nx)
    P[0] = 0.01
    
    for n in range(nt):  #iterate through time
    #def animate(j):
        Pn = P.copy()
        Tn = T.copy() ##copy the existing values of T into Tn
        for i in range(1, nx - 1):
       
            #T[i] = Tn[i] - lamb * ((dt[i]-dt[i-1]) /(dx[i]-dx[i-1])) * (Tn[i] - Tn[i-1]) +  ((dt[i]-dt[i-1]) / (dx[i]-dx[i-1])**2) * (1/Pe) * (Tn[i+1] - 2 * Tn[i] + Tn[i-1])
            T[i] = Tn[i] - lamb * (1/Length) * (dt /dx) * (Tn[i] - Tn[i-1]) + (1/Length**2)*(dt / dx**2) * (1/Pe) * (Tn[i+1] - 2 * Tn[i] + Tn[i-1])
            if T[i] >= 1.0 and T[i-1]<=1:
                constant = (f*rho_l*L*Q*Length)/(dx*i*kappa*rho_p*C_p)
                T[i] = Tn[i] - lamb * (1/Length) * (dt /dx) * (Tn[i] - Tn[i-1]) + (1/Length**2)*(dt / dx**2) * (1/Pe2) * ((2*Tn[i]-Tn[i-1]-constant*dx) - 2 * Tn[i] + Tn[i-1])
            
            theta[i] = (T[i] - 0)/(1-0)
            
            P[i] = ((kappa/(phi*Length**2))*(((1/Length**2)*((Pn[i]/Tn[i])-(Pn[i-1]/Tn[i-1]))*(Pn[i]-Pn[i-1]))/dx**2)+((1/Length**2)*(Pn[i]/Tn[i])*((Pn[i+1]-2*Pn[i]+Pn[i-1])/dx**2))*dt+(Pn[i]/Tn[i]))*T[i]
    
    #        line.set_data(numpy.linspace(0, length, nx), theta)
    #        
    #    return line,
        #pyplot.figure(1)
        #pyplot.semilogx(numpy.linspace(0, 2, nx), theta)
    pyplot.figure(2)
    pyplot.semilogx(numpy.linspace(0.001, length, nx), T)
    pyplot.semilogx(numpy.linspace(0.001, length, nx), P)
        #pyplot.plot(numpy.linspace(0, 2, nx), theta);
        
    # call the animator.  blit=True means only re-draw the parts that have changed.
    #anim = animation.FuncAnimation(fig, animate, init_func=init,
    #                               frames=20000, interval=20, blit=True)