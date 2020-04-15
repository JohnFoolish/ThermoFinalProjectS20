# -*- coding: utf-8 -*-
"""
Created on Tue Apr 14 17:37:27 2020

@author: John Lewis Corker and Somil Joshi
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import math
import pandas 
import argparse as arg


def main():
    
    stat_mech_model()
    
    
    
    
    
    
def stat_mech_model():
    """
    This function plots the results of the model calculated from statistical 
    mechanics. The model itself was derived by the Stat Mech textbook/Chesson
    Sipling.
    
    Inputs:
    -------
    None
    
    Returns:
    --------
    Graph: A matplotlib graph of S, T, state-space.    
    
    """

    #Fill in some dummy values for the Entropy Equation
    m = 1.0
    h = 6.626070e-34
    C = (np.e*m**(3/2)/h**3) * ((4 * np.pi * np.e)/3)**(3/2)
    v = 1.0
    N = 3.0
    n = 10.0e20

    #Populate list array with a range of different entropy values for the 
    #monatomic ideal gas
    S = []
    T = []
    V = []
    
    #Now calculate the entropy for the Monotomic ideal gas as it goes towards 
    #increasing temperature values and increasing volumes 
    for v in range(1, 110, 10):
        for t in range(50, 550, 50):
            
            #Chesson's big brain method
            #E_t = 1/2 * 1.380649e-23 * t
            #w = C*(v/N)*(E_t)**(3*N/2)
            
            #Dr. Rocklin's Google Method
            w = (v**N)/np.math.factorial(N) * ((2 * np.pi*m*1.380649e-23*t)/h**2)**(3*N/2)

            S.append(1.380649e-23* n* np.log(w))
            T.append(t)
            V.append(v)
    
    
    #plotting stuff
    
    S_arr = np.array(S)
    T_arr = np.array(T)
    V_arr = np.array(V)
    
    #fig, (ax1, ax2) = plt.subplots(1, 2, sharex = True)
    #color code based on stability
    
    fig = plt.figure()
    ax3 = fig.add_subplot(111, projection='3d')
    ax3.set_xlabel("Temperature [K]")
    ax3.set_ylabel("Volume")
    ax3.set_zlabel("Entropy")
    ax3.scatter(T_arr, V_arr, S_arr)
    
    
    #ax1.set_title("S-T State Space")
    #ax1.set_xlabel("Temperature [K]")
    #ax1.set_ylabel("Entropy")
    #ax1.plot(T_arr, S_arr, label = "Entropy for ST statespace")
    #ax1.legend()

    

    #ax2.set_title("V-T State Space")
    #ax2.set_xlabel("Temperature [K]")
    #ax2.set_ylabel("Volume")
    #ax2.plot(T_arr, V_arr, label = "Line of Isoschors")
    #ax2.legend()

    
main()