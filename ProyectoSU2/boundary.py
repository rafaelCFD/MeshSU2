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
import numpy as np

#N = 31
#Fx = 600

def boundaryNodos(Ny, N):
    """
    Boundary
    Ny = Number of coordinates
    N = Number of radial nodes
    """
    
    bounSU = np.zeros([N-1,6]) 
    cont = 0
    
    for i in range(0,N-1):
        bounSU[i,0] = 3
        bounSU[i,1] = i
        bounSU[i,2] = i + 1
        bounSU[i , 3] = 3
        bounSU[i , 4] = i + Ny - (N - 1)
        bounSU[i , 5] = i + Ny - N + 2
        cont += 1
        if cont == N-1:
            bounSU[i,1] = i
            bounSU[i,2] = i - N + 2
            bounSU[i ,4] = i + Ny - (N - 1)
            bounSU[i ,5] = Ny - i - 1
            
    cont = 0

    return bounSU
