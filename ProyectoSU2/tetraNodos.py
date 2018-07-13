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
#Ny = 600
def tetraNodos(Ny, N):
    """
    structured mesh
    Ny = Number of coordinates
    N = Number of radial nodes
    """
    nodSU = np.zeros([Ny-N+1, 6])   
    cont = 0
    for i in range(0, Ny-N+1):
        nodSU[i,0] = 9
        nodSU[i,1] = i
        nodSU[i,2] = i + 1
        nodSU[i,3] = i + N
        nodSU[i,4] = i + N - 1
        cont += 1
        if cont == N - 1:
          nodSU[i,1] = i
          nodSU[i,2] = i - N + 2
          nodSU[i,3] = i + 1
          nodSU[i,4] = i + N - 1
          cont = 0
        nodSU[i,5] = i # Columna para enumerar las celdas de la malla
    return nodSU
