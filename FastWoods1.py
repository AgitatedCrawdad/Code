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
#########################################################
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
    nt = 10000    #the number of timesteps we want to calculate
    Pe = 10.0
    Pe2 = 0.1
    lamb = 1.37
    #dt = sigma * dx**2 / nu #dt is defined using sigma ... more later!
    dt = 0.0001
    rho_l = 1000.0 #Density of liquid
    L = 2260.0     #Latent Heat of vaporization
    rho_p = 800.0
    C_p = 3000.0
    f = 0.2
    Length = 100.0
    kappa = 0.143E-6
    Q = (0.143*0.006)*10**-6
    #dx = numpy.logspace(-3,2,nx)
    #dt = numpy.logspace(-7,1,nx)
    #L = numpy.linspace(0,100,nx)
    T = numpy.ones(nx)*2      #a numpy array with nx elements all equal to 1.
    theta = numpy.ones(nx)*2
    theta[0] = 0
    #u[int(0.5 / dx):int(1.0 / dx + 1)] = 2  #setting u = 2 between 0.5 and 1 as per our I.C.s
    T[0] = 0.01
    
    Tn = numpy.ones(nx) #our placeholder array, Tn, to advance the solution in time
    
    numdx = numpy.arange(nx)
    
    for i in range(nt):  #iterate through time
    #def animate(j):
        #T[0:int(0.1*j)]=1
        Tn = T.copy() ##copy the existing values of T into Tn
        #for i in range(1, nx - 1):
            #T[i] = Tn[i] - lamb * ((dt[i]-dt[i-1]) /(dx[i]-dx[i-1])) * (Tn[i] - Tn[i-1]) +  ((dt[i]-dt[i-1]) / (dx[i]-dx[i-1])**2) * (1/Pe) * (Tn[i+1] - 2 * Tn[i] + Tn[i-1])
        T[1:-1] = Tn[1:-1] - lamb * (1/Length) * (dt /dx) * (Tn[1:-1] - Tn[0:-2]) + (1/Length**2)*(dt / dx**2) * (1/Pe) * (Tn[2:] - 2 * Tn[1:-1] + Tn[0:-2])
        theta = (T - 0)/(1-0)
        index = numpy.argmax(T>=1)
        if index > 1:
            index = index - 1
        constant = (f*rho_l*L*Q*Length)/(dx*kappa*rho_p*C_p)
#        T[index:-1] = Tn[index:-1] - lamb * (1/Length) * (dt /dx) * (Tn[index:-1] - Tn[index-1:-2]) + (1/Length**2)*(dt / dx**2) * (1/Pe2) * ((2 * Tn[index:-1]-Tn[index-1:-2]-constant*(1/numdx[index:-1])*dx) - 2 * Tn[index:-1] + Tn[index-1:-2])
#        T[index:index+1] = Tn[index:index+1] - lamb * (1/Length) * (dt /dx) * (Tn[index-1:index] - Tn[index-1:index]) + (1/Length**2)*(dt / dx**2) * (1/Pe2) * ((2 * Tn[index:index+1]-Tn[index-1:index]-constant*(1/numdx[index:index+1])*dx) - 2 * Tn[index:index+1] + Tn[index-1:index])
        T[index:index+1] = T[index:index+1] - lamb * (1/Length) * (dt /dx) * (T[index-1:index] - T[index-1:index]) + (1/Length**2)*(dt / dx**2) * (1/Pe2) * ((2 * T[index:index+1]-T[index-1:index]-constant*(1/numdx[index:index+1])*dx) - 2 * T[index:index+1] + T[index-1:index])

#        if T[T>0.65]>0.65:
#            constant = (f*rho_l*L*Q*Length)/(kappa*rho_p*C_p)
#            T[1:-1] = Tn[1:-1] - lamb * (1/L) * (dt /dx) * (Tn[1:-1] - Tn[0:-2]) + (1/L**2)*(dt / dx**2) * (1/Pe) * ((2*Tn[1:-1]-Tn[0:-2]-constant) - 2 * Tn[1:-1] + Tn[0:-2])
            
    #        line.set_data(numpy.linspace(0, length, nx), theta)
    #        
    #    return line,
        #pyplot.figure(1)
        #pyplot.semilogx(numpy.linspace(0, 2, nx), theta)

    pyplot.figure(2)
    pyplot.semilogx(numpy.linspace(0.001, length, nx), T)
    pyplot.axis([0.001, 100, 0,2 ])
    pyplot
        #pyplot.plot(numpy.linspace(0, 2, nx), theta);
        
    # call the animator.  blit=True means only re-draw the parts that have changed.
    #anim = animation.FuncAnimation(fig, animate, init_func=init,
    #                               frames=20000, interval=20, blit=True)