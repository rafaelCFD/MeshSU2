#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#########################################################################################################
#  \file mainSU2.py                                                                                     #
#  \Python script for creating a .su2 structured mesh of a 4-digit NACA airfoil with circular domain.   #
#  \author Rafael Medina Nogueron                                                                       #
#  \Instituto Politecnico Nacional                                                                      #
#  \rafaelmn@ipn.mx                                                                                     #
#  \version 1.0.                                                                                        #
#########################################################################################################
#import matplotlib.pyplot as plt
import numpy as np

def airfoil(c,N,m1,p1,t1):
    """
   Generation of the airfoil  NACA 4D
   c     = Chord line
   N     =  number of points (debe ser inpar)
   m1    = x is the position along the chord 
   p1    = y_{t} is the half thickness at a given value of x
   t1    = t is the maximum thickness as a fraction of the chord (so 100 t gives the last two digits in the NACA 4-digit denomination)
    
   x,y coordenadas cartesianas
   """ 
    #c  = 2.
    #N  = 21 
    #m1 = 2. 
    #p1 = 5. 
    #t1 = 12.
    
    m = m1/100
    p = p1/10
    t = t1/100
    
    a0 =  0.2969 
    a1 = -0.1260
    a2 = -0.3516 
    a3 =  0.2843 
    a4 = -0.1036 # Bs abierto a4=-0.1015;
    ele  = int((N+1)/2)
    xc = np.linspace(0, 1, ele)
    yc = np.ones(ele)
    yc_x = np.ones(ele)
    theta = np.ones(ele)
    
    # Camber line
    for i in range(ele):
        if xc[i] >= 0 and xc[i] < p:
            yc[i] = m/p**2*(2*p*xc[i]-xc[i]**2)
            yc_x[i] = ((2*m)/(p**2))*(p-xc[i])
        elif xc[i] >= p and xc[i] <= 1:
            yc[i] = m/(1-p)**2*(1-2*p+2*p*xc[i]-xc[i]**2)
            yc_x[i] = ((2*m)/((1-p)**2))*(p-xc[i])
        theta[i] = np.arctan(yc_x[i])
    
    # Thickness distribution
    yt = 5*t*((a0*xc**(1/2)) + (a1*xc) + (a2*xc**2) + (a3*xc**3) + (a4*xc**4))
    
    # Extrados nodes
    xu = xc[:] - yt[:]*np.sin(theta)
    yu = yc[:] + yt[:]*np.cos(theta)
    
    # Intrados nodes
    xl = xc[:] + yt[:]*np.sin(theta)
    yl = yc[:] - yt[:]*np.cos(theta)
    

    cordxy = np.zeros([2*ele-1,2])
    for j in range(ele-1,0,-1):
        i = ele - j - 1
        cordxy[i,0] = xu[j]
        cordxy[i,1] = yu[j]
    for j in range(0,ele):
        i = ele + j - 1
        cordxy[i,0] = xl[j]
        cordxy[i,1] = yl[j]

    cordxy = c * cordxy
    cordxy[:,0] = cordxy[:,0] - c * 0.5
    
    # plot airfoil
    #plt.plot(cordxy[:,0],cordxy[:,1])
    #plt.axis('equal')
    #plt.grid()
    
    np.savetxt('perfilCoordenadas.txt', cordxy)
    
    return cordxy

#airfoil(c,N,m1,p1,t1)
#airfoil(2.,21,2.,5.,12.)
