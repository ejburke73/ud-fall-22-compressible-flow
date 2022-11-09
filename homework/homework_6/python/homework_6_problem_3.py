# Compressible Flow
# AEE 553
# Homework 6 - Problem 2
# Evan Burke

import numpy as np
from matplotlib import pyplot as plt
import shocks as ns
import oblique as os
import isentropic as isen
from scipy.optimize import fsolve

pt = 20.408 * 10**6 # MPa to Pa
T_i = 3600 # static temp at inlet of CD nozzle, K
A_i = 0.21 # area at inlet of CD nozzle, m^2
A_th = 0.054 # area at throat of CD nozzle, m^2
A_e = 4.17 # area at exit of CD nozzle, m^2
R = 287
gamma = 1.2 # this is different!!!

def mass_flow(pt=None, A_star=None, Tt=None, gamma=1.4, R=287):
    mdot = pt*A_star / Tt**0.5 * (gamma/R * (2/gamma+1)**((gamma+1)/(gamma-1)))**0.5
    print(f'Choked Mass Flow Rate = {mdot} kg/s')

def A_from_A_star(A_star=None,M=None,gamma=1.4):
    A = ((A_star**2/M**2) * (2/(gamma+1) * (1 + (gamma-1)/2 * M**2 ))**((gamma+1)/(gamma-1)))**0.5
    print(f'Area for M = {M}: {A} m^2')
    return A  

altitudes = np.linspace(0,20000,num=20001, endpoint=True)
print(altitudes)

pressures = [101325*(1-(2.25577*10**(-5)*h))**5.25588 for h in altitudes]
print(pressures[0:10])

def SSME_thrust(mdot=None,u_e=None,p_e=None,p_amb=None,A_e=None):
    thrust = mdot * u_e + (p_e-p_amb)*A_e
    pass
# To get mdot:
# Need pt, A_th, Tt, gamma, R
# pt is given, Tt unknown, M_th = 1, A_i/A_th known
# Solve for M_i from A_i/A_th given that M_th = 1

def A_A_star(M=None,A_A_star=None,gamma=1.4):
    eq = ((1/M**2) * (2/(gamma+1) * (1 + (gamma-1)/2 * M**2 ))**((gamma+1)/(gamma-1)))**0.5 - A_A_star
    return eq

M_i = float(fsolve(A_A_star,x0=0.01,args=(A_i/A_th,gamma)))
print(f'Nozzle Inlet Mach = {M_i}')
Tt = isen.get_total_temperature(M=M_i,T=T_i,gamma=gamma)
print(f'Total temperature at nozzle inlet = {Tt} K')

mdot = mass_flow(pt=pt,A_star=A_th,Tt=Tt,gamma=gamma,R=R)

# mdot, A_e, p_ambs known, need u_e and p_e
# Use A/A* relationship to get exit Mach number, exit sonic velocity, exit velocity
