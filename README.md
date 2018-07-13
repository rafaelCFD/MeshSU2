# MeshSU2
Create a mesh of a airfoil

#########################################################################################################
#  \file mainSU2.py                                                                                     #
#  \Python script for creating a .su2 structured mesh of a 4-digit NACA airfoil with circular domain.   #
#  \author Rafael Medina Nogueron                                                                       #
#  \Instituto Politecnico Nacional                                                                      #
#  \rafaelmn@ipn.mx                                                                                     #
#  \version 1.0.                                                                                        #
#########################################################################################################

Run the file mainSU2 in python 3

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
