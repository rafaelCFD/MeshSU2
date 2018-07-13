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
import matplotlib.pyplot as plt
import numpy as np
#import math

import airfoil as af
import mesh as mh
import tetraNodos as tn
import boundary as bd

#########################################################################################################
########################################### Input data ##################################################
#########################################################################################################
N  = 101        # Airfoil number of points
M  = 50        # 20
c  = 1.        # Chord [m]
m  = 2.        # Is the maximum camber (100 m is the first of the four digits),
p  = 5.        # Is the location of maximum camber (10 p is the second digit in the NACA xxxx description).
t  = 12.       # Is the maximum thickness as a fraction of the chord 
                # (so t gives the last two digits in the NACA 4-digit denomination divided by 100)
#########################################################################################################
######################################################################################################### 
#########################################################################################################

af.airfoil( c , N , m , p , t )

cordxy = np.loadtxt('perfilCoordenadas.txt',float)

[X, Y] = mh.mesh( 1. , M , cordxy )


########################### Plot mesh & airfoil ###############################
plt.plot(X,Y,"b",X.T,Y.T,"g")
plt.title('Airfoil NACA %d%d%d' % (m, p, t))
plt.xlabel('x')
plt.ylabel('y')
plt.grid()
plt.axis("equal")  

##################### Coordinates, Mesh & Boundary ############################
X = X[:-1,:]    # remove the last element from the column
X = X.T.ravel() # Convert an array in a column
Y = Y[:-1,:]
Y = Y.T.ravel() 

cordSU      = np.zeros([Y.size,3])
cordSU[:,0] = X; cordSU[:,1] = Y
for i in range(0,Y.size):
    cordSU[i,2] = i

nodSU = tn.tetraNodos(Y.size, N)

bounSU = bd.boundaryNodos(Y.size, N)

##################### Generate the file .SU2 ##################################
n     = 2                                  # Dimenciones
nFron = 2                                  # Numero de fronteras

f = open ("mesh_NACA0012_inv.su2", "w")

f.write("NDIME= %d \n" % (n))
f.write("NELEM= %d \n" % (nodSU.shape[0]))
for i in range(0,nodSU.shape[0]):
    f.write("%d \t %d \t %d \t %d \t %d \t %d \n" % (nodSU[i,0],nodSU[i,1],nodSU[i,2],nodSU[i,3],nodSU[i,4],nodSU[i,5]))

f.write("NPOIN= %d \n" % (cordSU.shape[0]))
for i in range(0,cordSU.shape[0]):
    f.write("%f \t %f \t %d \n" % (cordSU[i,0],cordSU[i,1],cordSU[i,2]))

f.write("NMARK= %d \n" % (nFron))
f.write("MARKER_TAG= airfoil \n")
f.write("MARKER_ELEM= %d \n" % (bounSU.shape[0]))
for i in range(0,bounSU.shape[0]):
    f.write("%d \t %d \t %d \n" % (bounSU[i,3],bounSU[i,4],bounSU[i,5]))

f.write("MARKER_TAG= farfield \n")
f.write("MARKER_ELEM= %d \n" % (bounSU.shape[0]))
for i in range(0,bounSU.shape[0]):
    f.write("%d \t %d \t %d \n" % (bounSU[i,0],bounSU[i,1],bounSU[i,2]))
f.close()
