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
from graph_class import graph_data


def main():
    
    Ideal_statmech_model()
    Ideal_classical_model()
    VDW_classical_model()
    ##Add command line stuff so that you can adjust the axis to observe 
    #different kind of changing variables (ie: E-V state space, V-P, etc)
    
    
def VDW_classical_model():
    """
    This function plots the results of the VDW model for entropy calculated 
    using classical Thermodynamics.

    Inputs:
    -------
    None

    Returns:
    --------
    Graph: A matplotlib graph of S, T, state-space.
    Equations:

    Ideal Gas:

    PV = NRT
    P = pressure
    V = volume
    N = moles of substance = N_particles / N_a
    R = ideal gas constant = N_a * K_b
    T = temperature in Kelvin

    U = cNRT
    U = energy
    c = (# degrees of freedom / 2)

    Molar Volume: v = V / N
    Molar Energy: u = U / N
    Molar Entropy: s = S / N

    Entropy is extensive:
    S(U, V, N) = N * S(U / N, V / N, 1) = N * s(u, v)

    ds = (1 / T) * du + (P /T) dv
    s = cR*log(u) + R*log(v) + constant

    VDW Gas:
    P = (RT) / (v -b) - a / (v^2)
    R = ideal gas constant = N_a * K_b
    a, b = material specific constants
    v = molar volume
    ...
    Molar Entropy: s(u, v) = cR*log(u + a/v) +R*log(v - b) + constant
    Entropy: S = s * N
    """

    # constants
    c = 2 / 2
    R = 8.3145
    # Van der Waals constants for hydrogen gas
    a = 0.2453
    b = 0.02651
    N = 1.0

    # Populate list array with a range of different entropy values
    Slist = []
    Ulist = []
    Vlist = []
    Tlist = []
    
    graph_data_class = graph_data("Van Der Waals Gas", "Volume [L]", "Temperature [K]", "Entropy [J/K]")
    
    
    

    # Now calculate the entropy for VDW gas as it goes towards
    # increasing temperature values
    for V in range(1, 101, 5):
        v = V / N
        for T in range(273, 374, 5):
            U = (3 / 2) * R * T
            u = U /N
            s = (c*R)*np.log(u + a/v) + R*np.log(v - b)
            
            #Add the data to the list
            Slist.append(s * N)
            Ulist.append(U)
            Vlist.append(V)
            Tlist.append(T)
        
            
            if (len(Slist) < 2):
                stable = False
            else:
                
                #At constant Volume, stability is found through dE/dT > 0.
                C_v = (Ulist[1] - Ulist[0])/(Tlist[1] - Tlist[0])
                
                if(C_v > 0):
                    stable = True
                else:
                    stable = False
                
                #Reset the list so that we accurately calculate the stability 
                #for each point
                tmpU = Ulist[1]
                tmpT = Tlist[1]
                Ulist.clear()
                Tlist.clear()
                Ulist.append(tmpU)
                Tlist.append(tmpT)
            
            
            
            graph_data_class.plot_point(v, T, s, stable)
            
            
    graph_data_class.display()


def Ideal_classical_model():
    """
    This function plots the results of the Ideal gas model for entropy
    calculated using classical Thermodynamics.

    Inputs:
    -------
    None

    Returns:
    --------
    Graph: A matplotlib graph of S, T, state-space.
    Equations:

    Ideal Gas:

    PV = NRT
    P = pressure
    V = volume
    N = moles of substance = N_particles / N_a
    R = ideal gas constant = N_a * K_b
    T = temperature in Kelvin

    U = cNRT
    U = energy
    c = (# degrees of freedom / 2)

    Molar Volume: v = V / N
    Molar Energy: u = U / N
    Molar Entropy: s = S / N

    Entropy is extensive:
    S(U, V, N) = N * S(U / N, V / N, 1) = N * s(u, v)

    ds = (1 / T) * du + (P /T) dv
    s = cR*log(u) + R*log(v) + constant
    """

    # constants
    c = 2 / 2
    R = 8.3145
    N = 1.0


    
    graph_data_class = graph_data("Ideal Gas", "Volume [L]", "Temperature [K]", "Entropy [J/K]")

    # Populate list array with a range of different entropy values
    Slist = []
    Ulist = []
    Vlist = []
    Tlist = []

    # Now calculate the entropy for VDW gas as it goes towards
    # increasing temperature values
    for V in range(1, 101, 5):
        v = V / N
        for T in range(273, 374, 5):
            U = (3 / 2) * R * T
            u = U /N
            s = (c*R)*np.log(u) + R*np.log(v)
            Slist.append(s * N)
            Ulist.append(U)
            Vlist.append(V)
            Tlist.append(T)

            if (len(Slist) < 2):
                stable = True
            else:
                
                #At constant Volume, stability is found through dE/dT > 0.
                C_v = (Ulist[1] - Ulist[0])/(Tlist[1] - Tlist[0])
                
                if(C_v > 0):
                    stable = True
                else:
                    stable = False
                
                #Reset the list so that we accurately calculate the stability 
                #for each point
                tmpU = Ulist[1]
                tmpT = Tlist[1]
                Ulist.clear()
                Tlist.clear()
                Ulist.append(tmpU)
                Tlist.append(tmpT)
            
            
            
            graph_data_class.plot_point(v, T, s, stable)
            
            
    graph_data_class.display()

    
def Ideal_statmech_model():
    """
    This function plots the entropy of an ideal gas system from the Statistical
    Mechanics model
    
    Inputs:
    -------
    None
    
    Returns:
    --------
    Graph: A matplotlib graph of S, T, state-space.   
    
    Equations:
    ----------
    m = mass of particle
    h = Plank's constant
    C = ((e * m) ^ (3/2) ) / h^3 * ((4pi * e) / 3) ^ (3/2)   
    v = Volume (Liters)
    N = number of particles
    n = scaling factor
    t = temperature
    U = Total Thermal energy given by 1/2 * K_b * t
    
    w = number of microstates 
    
    
    
    """

    #Fill in some dummy values for the Entropy Equation
    m = 1.0
    h = 6.626070e-34
    v = 1.0
    N = 3.0
    n = 10.0e20
    
    
    graph_data_class = graph_data("Stat Mech Ideal Gas", "Volume [L]", "Temperature [K]", "Entropy [J/K]")

    #Populate list array with a range of different entropy values for the 
    #monatomic ideal gas
    Slist = []
    Tlist = []
    Vlist = []
    Ulist = []
    
    #Now calculate the entropy for the Monotomic ideal gas as it goes towards 
    #increasing temperature values and increasing volumes 
    for v in range(1, 101, 5):
        for t in range(273, 374, 5):
            
            #Calculate the Thermal Energy of the system
            U = 1/2 * 1.380649e-23 * t
            
            w_N1 = ((np.e * v)/(N*h**3))**N
            
            w_N2 = ((4*np.pi*np.e*m*U)/3*N)**(3*N/2)
            
            #Calculate the number of microstates for the system
            w_N = w_N1 * w_N2
            
            scale_factor = 1#n         
            
            w_N_corrected = w_N * scale_factor
            
            s = 1.380649e-23*np.log(w_N_corrected)
            
            #s = w_N_corrected/1.380649e-23
            Slist.append((s))
            Tlist.append(t)
            Vlist.append(v)
            Ulist.append(U)
    
            if (len(Slist) < 2):
                stable = True
            else:
                
                #At constant Volume, stability is found through dE/dT > 0.
                C_v = (Ulist[1] - Ulist[0])/(Tlist[1] - Tlist[0])
                
                if(C_v > 0):
                    stable = True
                else:
                    stable = False
                
                #Reset the list so that we accurately calculate the stability 
                #for each point
                tmpU = Ulist[1]
                tmpT = Tlist[1]
                Ulist.clear()
                Tlist.clear()
                Ulist.append(tmpU)
                Tlist.append(tmpT)
            
            
            
            graph_data_class.plot_point(v, t, s, stable)
            
            
    graph_data_class.display()
    

    
main()