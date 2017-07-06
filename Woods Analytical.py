# -*- coding: utf-8 -*-
"""
Created on Thu Jun 29 11:31:28 2017

@author: dfranken
"""
import numpy as np
from scipy import integrate
from matplotlib import pyplot as plt

#Integration Equation
def f(eta):
    return eta**(float(lamb)*float(beta)-1)*np.exp(-(eta**2))

def bounds_x(i):
    return [0, eta1[i]]
def bounds_x2(i):
    return [0, np.inf]

#Constants
r = 1
r2 = np.linspace(0.001,100,1000)
exponent = 1
k = 0.143*10**(-exponent)
t1 = np.linspace(0.1,1000,1000)
t2 = 10
Q = 0.5*10**(-exponent)
phi = 0.5
T0 = 60
T1 = 20
lamb = 1
beta = Q/k
T = np.zeros(len(t1))

omega = (Q/(2*(phi*k))**(0.5))

eta1 = r/(2*(float(k)*t1)**(0.5))
eta2 = r2/(2*(float(k)*t2)**(0.5))

#Equation 3 from Liquid and Vapor Flow in Superheated Rock Paper for cold liquid injected into hot liquid
def Ex1():
    for i in range(len(t1)):    
        integrand, residual = integrate.nquad(f,[bounds_x(i)])
        deltaT, residual2 = integrate.nquad(f,[bounds_x2(i)])
        T[i] = T0 + ((T1-T0)/(deltaT))*integrand
    plt.figure(1)
    plt.semilogx(eta1,T)
#    plt.figure(2)
#    plt.plot(eta,T)


def Ex2():
    for i in range(len(t1)):    
        integrand, residual = integrate.nquad(f,[bounds_x(i)])
        deltaT, residual2 = integrate.nquad(f,[bounds_x2(i)])
        T[i] = T0 + ((T1-T0)/(deltaT))*integrand
    plt.figure(1)
    plt.semilogx(eta2,T)