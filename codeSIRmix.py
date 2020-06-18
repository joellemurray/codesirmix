#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 5 14:48:58 2020

@author: jmurray
"""
# uses the euler method to numerically solve coupled first order ODE's
# for the SIR ODE model using the euler method
# S = Susceptible, I = Infected, R = Removed (recovered/deceased/immune)
# dS/dt = -beta*S*I/N, dI/dt = beta*S*I-gamma*I, dR/dt=gamma*I
# where S+I+R=constant=total population
# https://www.davidketcheson.info/2020/03/19/SIR_Estimating_parameters.html
# https://www.davidketcheson.info/2020/03/19/SIR_predictions.html
# basic assumption beta = 0.25 and gamma = 0.06
# added in 2 different susceptible populations: S1 and S2
# S1 uses physical distancing to change the beta to q*beta where q = 0 for 
# full PD (total isolation) and q=1 for no PD
# while S2 never uses PD

# importing packages
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy as sp
import math

#initializing parameters
# time step (day)
dt = 0.01
# average number of people that come within infection range of an infected individual
# per day
beta = 0.25
# remove probability per day = 1/(recovery time)
# approximately 16 days (varies from 14 days to a month)
gamma = 0.05

# total population
N = 1E6
# total time (days)
tottime = 365
# initial percent removed (immune)
pr = 0.0
# initial percent infected
pi = 0.01
# initial percent susceptible
ps = 1-pr-pi
# percent isolators (physical distance = PD)
ppd = 0.5
# plotting
fig = plt.figure()
for i in range(0,11):
# social distance factor (0 = total isolation, 1 = no isolation)
    q = 0.1*(i-1)
    
    niter = int(math.ceil(tottime/dt))
    t = np.arange(0, tottime, dt)   
    S = np.zeros(niter)
    S1 = np.zeros(niter)
    S2 = np.zeros(niter)
    I = np.zeros(niter)

# initial populations
 
    S[0] = ps*N
    S1[0]= ppd*S[0]
    S2[0]= S[0]-S1[0]
    I[0] = pi*N

    
    for j in range(niter-1):        
        dS1dt=-q*beta/N*S1[j]*I[j]
        dS2dt=-beta/N*S2[j]*I[j]
        dIdt=-dS1dt-dS2dt-gamma*I[j]
        S1[j+1] = S1[j] + dt*dS1dt 
        S2[j+1] = S2[j] + dt*dS2dt
        I[j+1] = I[j] + dt*dIdt  
        
    S=S1+S2   
    R = N-S-I

#plt.plot(t, S, 'k', label = 'susceptible')
    hfont = {'fontname':'DejaVu Sans'}
    mpl.rcParams["font.size"] = 12
    
    plt.plot(t, I, label = 'infected')
#plt.plot(t, R, 'b', label = 'recovered')
#plt.plot(t, R+S+I, 'y', label = 'total')
#    plt.gca().legend(('susceptible','infected','recovered','total'))
    plt.title(r'codesirmix model with N = 1E6, $\beta$  = 0.25, $\gamma$ = 0.05')
    plt.xlabel('days elapsed since 1 percent of the population became infected')
    plt.ylabel('total infected population')



plt.annotate('q = 0.0', ha = 'center', va = 'bottom', xytext = (1.3E2, 1.1E5), 
             xy = (60, 0.8E5),arrowprops = {'facecolor' : 'black'})
plt.annotate('q = 1.0', ha = 'center', va = 'bottom', xytext = (1.2E2, 4.0E5), 
             xy = (40, 4.6E5),arrowprops = {'facecolor' : 'black'})
plt.annotate('50% PD population', xy=(0.6, 0.8), xycoords='axes fraction', **hfont)
plt.annotate('q = 0.0 -> full PD', xy=(0.6, 0.74), xycoords='axes fraction')
plt.annotate('q = 1.0 -> no PD', xy=(0.6, 0.68), xycoords='axes fraction')
plt.annotate(r'$\frac{\partial S_1}{\partial t} = - \frac{q \beta}{N} S_1 I$', 
             xy=(0.6, 0.58), xycoords='axes fraction')
plt.annotate(r'$\frac{\partial S_2}{\partial t} = -\frac{\beta}{N} S_2 I$', 
             xy=(0.6, 0.48), xycoords='axes fraction')
plt.annotate(r'$S = S_1 + S_2$', xy=(0.6, 0.38), xycoords='axes fraction')
plt.show()

