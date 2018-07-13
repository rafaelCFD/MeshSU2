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

def mesh(c,M,cordxy):
    """
    c     = Chord line
    M     =  number of points (radial eta)
    Mesh generation (Jacobi Method)
    xi, eta  = generalized coordinates
    X,Y      = coordinates x,y
    it       = iterations
    it_max   = maximum iterations
    d_dx     = Variación de la coordenada x con respecto a la iteración anterior
    d_dy     = Variación de la coordenada y con respecto a la iteración anterior
    err      = Error
    Q        = Función para generar atracción en una línea coordenada
    ek,Ak,Ck = Parámetros que regulan la atracción de la línea coordenada
    w        = Factor de relajación, controla la velocidad de convergencia de la solución 
    External boundary (circle)
    r     = radio 
    N     = discretization circle
    """
    
    # Internal boundary
    
    #c = 2.
    #M = 30
    N = cordxy.shape[0] # airfol coordinates (airfoil.py)
    X = np.zeros([N,M])
    Y = np.zeros([N,M])
    
    X[:,M-1] = cordxy[:,0]
    Y[:,M-1] = cordxy[:,1]    
    
    theta = np.linspace(0,2*np.pi,N)
    r = c
    
    # External boundary
    R = 5 * r
    
    X[:,0] = R * np.cos(theta) 
    Y[:,0] = R * np.sin(theta)
    
    # Spacing between axis eta
    for j in range(1,M-1):
        X[:,j] = X[:,j-1] + (X[:,M-1]-X[:,0])/(M-1)
        Y[:,j] = Y[:,j-1] + (Y[:,M-1]-Y[:,0])/(M-1)
        
    # Loop Poisson funtion
    ek = M; Ak = 1; Ck = 0.2
    eta = np.arange(0,M) 
    xi = np.arange(0,N)
    
    it = 0 
    it_max = 1000
    err    = 1.e-8 
    d_dx   = 1
    d_dy   = 1   
    A = np.zeros([N,M])
    B = np.zeros([N,M])
    C = np.zeros([N,M])
    I = np.zeros([N,M])
    Q = np.zeros([N,M])
    P = np.zeros([N,M])    
    
    while d_dx > err and d_dy > err and it < it_max:
        it = it + 1
        Xn = X.copy()
        Yn = Y.copy()
        
        # Increased mesh density on the surface
        for i in range(0,N):
            for j in range(0,M):
                Q[i,j] = Ak*np.exp(-Ck*abs(eta[j]-ek))
                P[i,j] = 0
                
        # Internal nodes of the mesh (Poisson equation)
        for j in range(1,M-1):
            for i in range(1,N):
                if i == N-1:
                    A[N-1,j] = 0.25 * ((X[i,j+1]-X[i,j-1])**2+(Y[i,j+1]-Y[i,j-1])**2)
                    B[N-1,j] = (0.25 * ((X[1,j]-X[i-1,j])*(X[i,j+1]-X[i,j-1])) +
                               0.25*((Y[1,j]-Y[i-1,j])*(Y[i,j+1]-Y[i,j-1])))
                    C[N-1,j] = 0.25*(X[1,j]-X[i-1,j])**2+0.25*(Y[1,j]-Y[i-1,j])**2
                    I[N-1,j] = 0.25*((X[1,j]-X[i-1,j])*(Y[i,j+1]-Y[i,j-1])+ 
                               (Y[1,j]-Y[i-1,j])*(X[i,j+1]-X[i,j-1]))
        
                    X[N-1,j] = (0.5/(A[i,j]+C[i,j]))*(A[i,j]*(X[1,j]+X[i-1,j]) 
                               +C[i,j]*(X[i,j+1]+X[i,j-1])-0.5*B[i,j]*(X[1,j] 
                               -X[1,j-1]-X[i-1,j+1]+X[i-1,j-1])+0.5*I[i,j]**2*(P[i,j]*(X[1,j] 
                               -X[i-1,j])+Q[i,j]*(X[i,j+1]-X[i,j-1])))
                    Y[i,j]   = (0.5/(A[i,j]+C[i,j]))*(A[i,j]*(Y[1,j]+Y[i-1,j]) 
                               +C[i,j]*(Y[i,j+1]+Y[i,j-1])-0.5*B[i,j]*(Y[1,j+1] 
                               -Y[1,j-1]-Y[i-1,j+1]+Y[i-1,j-1])+0.5*I[i,j]**2*(P[i,j]*(Y[1,j] 
                               -Y[i-1,j])+Q[i,j]*(Y[i,j+1]-Y[i,j-1])))
                else:
                    A[i,j]   = 0.25*((X[i,j+1]-X[i,j-1])**2+(Y[i,j+1]-Y[i,j-1])**2)
                    B[i,j]   = (0.25*((X[i+1,j]-X[i-1,j])*(X[i,j+1]-X[i,j-1])) 
                               +0.25*((Y[i+1,j]-Y[i-1,j])*(Y[i,j+1]-Y[i,j-1])))
                    C[i,j]   = 0.25*(X[i+1,j]-X[i-1,j])**2+0.25*(Y[i+1,j]-Y[i-1,j])**2
                    I[N-1,j] = 0.25*((X[i+1,j]-X[i-1,j])*(Y[i,j+1]-Y[i,j-1]) 
                               +(Y[i+1,j]-Y[i-1,j])*(X[i,j+1]-X[i,j-1]))
        
                    X[i,j]   = (0.5/(A[i,j]+C[i,j]))*(A[i,j]*(X[i+1,j]+X[i-1,j]) 
                               +C[i,j]*(X[i,j+1]+X[i,j-1])-0.5*B[i,j]*(X[i+1,j+1] 
                               -X[i+1,j-1]-X[i-1,j+1]+X[i-1,j-1])+0.5*I[i,j]**2 *(P[i,j]*(X[i+1,j] 
                               -X[i-1,j])+Q[i,j]*(X[i,j+1]-X[i,j-1])))
                    Y[i,j]   = (0.5/(A[i,j]+C[i,j]))*(A[i,j]*(Y[i+1,j]+Y[i-1,j]) 
                               +C[i,j]*(Y[i,j+1]+Y[i,j-1])-0.5*B[i,j]*(Y[i+1,j+1] 
                               -Y[i+1,j-1]-Y[i-1,j+1]+Y[i-1,j-1])+0.5*I[i,j]**2*(P[i,j]*(Y[i+1,j] 
                               -Y[i-1,j])+Q[i,j]*(Y[i,j+1]-Y[i,j-1])))
                
                   
        X[0,:] = X[N-1,:]
        Y[0,:] = Y[N-1,:]
        d_dx = np.max(np.abs(X-Xn))
        d_dy = np.max(np.abs(Y-Yn))
    # plot Mesh
    #plt.plot(X,Y,"b",X.T,Y.T,"g")
    #plt.grid()
    #plt.axis("equal")      
    # otra forma de graficar
    #for j in range(0,M):
    #    for i in range(0,N):
    #        plt.plot(X[i,:],Y[i,:])#; hold on;axis equal;
    #        plt.plot(X[:,j],Y[:,j])#; hold on;axis equal;
    #plt.grid()
    #plt.axis("equal")

    return X,Y

#cordxy = np.loadtxt('perfilCoordenadas',float)
#mesh(1.,50,cordxy)