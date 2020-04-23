import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import math
import pandas
import argparse as arg


def main():
    VDW_classical_model()
    Ideal_classical_model()

def VDW_classical_model():
    """
    This function plots the results of the VDW model for entropy calculated using classical Thermodynamics.

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

    # Now calculate the entropy for VDW gas as it goes towards
    # increasing temperature values
    for V in range(1, 100, 1):
        v = V / N
        for T in range(273, 373, 10):
            U = (3 / 2) * R * T
            u = U /N
            s = (c*R)*np.log(u + a/v) + R*np.log(v - b)
            Slist.append(s * N)
            Ulist.append(U)
            Vlist.append(V)

    # plotting stuff

    S_arr = np.array(Slist)
    U_arr = np.array(Ulist)
    V_arr = np.array(Vlist)

    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111, projection='3d')
    ax1.set_xlabel("Energy [J]")
    ax1.set_ylabel("Volume [L]")
    ax1.set_zlabel("Entropy [J / K]")
    ax1.scatter(U_arr, V_arr, S_arr)
    plt.title("VDW Gas Entropy")

    plt.show()

def Ideal_classical_model():
    """
    This function plots the results of the Ideal gas model for entropy calculated using classical Thermodynamics.

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

    # Populate list array with a range of different entropy values
    Slist = []
    Ulist = []
    Vlist = []

    # Now calculate the entropy for VDW gas as it goes towards
    # increasing temperature values
    for V in range(1, 100, 1):
        v = V / N
        for T in range(273, 373, 10):
            U = (3 / 2) * R * T
            u = U /N
            s = (c*R)*np.log(u) + R*np.log(v)
            Slist.append(s * N)
            Ulist.append(U)
            Vlist.append(V)

    # plotting stuff

    S_arr = np.array(Slist)
    U_arr = np.array(Ulist)
    V_arr = np.array(Vlist)

    fig = plt.figure()
    ax3 = fig.add_subplot(111, projection='3d')
    ax3.set_xlabel("Energy [J]")
    ax3.set_ylabel("Volume [L]")
    ax3.set_zlabel("Entropy [J / K]")
    ax3.scatter(U_arr, V_arr, S_arr)
    plt.title("Ideal Gas Entropy")

    plt.show()




main()