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

pt = 4000*1000 # Pa
Tt = 500 # K
Me = 6 # Design exit Mach
D_th = 4.114 # throat diameter, inviscid, inches
D_th_r = 3.71 # throat diameter, real, viscous, inches
R = 287
gamma = 1.4

# Convert diameters to meters
D_th = D_th * 0.0254
D_th_r = D_th_r * 0.0254

print(f'Throat Diameter (Inviscid) = {D_th}')
print(f'Throat Diameter (Viscous) = {D_th_r}')

A_th = np.pi * D_th**2/4
print(f'A* = {A_th}')

def mass_flow(pt=None, A_star=None, Tt=None, gamma=1.4, R=287):
    mdot = pt*A_star / Tt**0.5 * (gamma/R * (2/(gamma+1))**((gamma+1)/(gamma-1)))**0.5
    print(f'Choked Mass Flow Rate = {mdot} kg/s')
    return mdot

def A_from_A_star(A_star=None,M=None,gamma=1.4):
    A = ((A_star**2/M**2) * (2/(gamma+1) * (1 + (gamma-1)/2 * M**2 ))**((gamma+1)/(gamma-1)))**0.5
    print(f'Area for M = {M}: {A} m^2')
    return A  
    
# Part A
mdot = mass_flow(pt=pt, A_star=A_th, Tt=Tt, gamma=1.4, R=287)

# Part B
A_exit = A_from_A_star(A_star=A_th,M=Me,gamma=1.4)

# Part C
p_exit = isen.get_static_pressure(M=Me,p_t=pt)
T_exit = isen.get_static_temperature(M=Me,T_t=Tt)

# Part D
def A_A_star(M=None,A_A_star=None,gamma=1.4):
    eq = ((1/M**2) * (2/(gamma+1) * (1 + (gamma-1)/2 * M**2 ))**((gamma+1)/(gamma-1)))**0.5 - A_A_star
    #print(f'A/A* = {A_A_star}')
    return eq

M_e_sub = float(fsolve(A_A_star,x0=0.01,args=(A_exit/A_th)))
print(f'Subsonic Exit Mach = {M_e_sub}')
p_e_sub = isen.get_static_pressure(M=M_e_sub,p_t=pt)

# Part E

# Normal shock will stand at nozzle exit when 
# static pressure is equal to the static pressure
# across a normal shock at the nozzle design condition
p_exit_ns = ns.get_static_pressure_normal_shock(M1=Me,p1=p_exit)
print(f'Back Pressure for which exit NS= {p_exit_ns}')

# Part F
print(f'Back Pressure below which no shocks in nozzle = {p_exit_ns}')

# Part G
# Range of back pressures for which there are oblique shocks
# in nozzle exhaust
# pe < pb
print(f'Range of back pressures for oblique shocks: {p_exit} < p_b < {p_exit_ns}')

# Part H
# Range of back pressures for expansion waves
# pb < pe

print(f'Range of back pressures for expansion waves: p_b < {p_exit}')

# Part I

A_avg = (A_th + A_exit)/2
print(f'Average nozzle area = {A_avg} m^2')
M_avg = float(fsolve(A_A_star,x0=1.5,args=(A_avg/A_th)))
print(f'M_avg = {M_avg}')
p_avg = isen.get_static_pressure(M=M_avg,p_t=pt)
p_avg_NS = ns.get_static_pressure_normal_shock(M1=M_avg,p1=p_avg)

# Part J
T_tank = 295 
tank_vol = 4000 # gal
tank_vol = tank_vol * 0.00378541 # m^3

ts = np.linspace(0,300,num=301,endpoint=True)
pbs = [mdot/tank_vol*t*R*T_tank for t in ts]

t_NS = p_exit_ns * tank_vol / (mdot * R * T_tank)
print(f'Time to fill receiving tanks to yield exit NS = {t_NS}')

# Part K
mass_used =mdot*t_NS
print(f'Mass used until NS = {mass_used}')
rho_tube = isen.get_static_density(p=pt,T=Tt)
D_driver = 9.75
D_driver = 0.24765
A_driver = np.pi*D_driver**2/4
volume_used = mass_used/rho_tube
L_used = mass_used/(rho_tube*A_driver)

print(f'Volume used = {volume_used}')
print(f'Driver Tube cross sectional area = {A_driver}')
print(f'Length used = {L_used} m = {L_used*3.28084} ft')
