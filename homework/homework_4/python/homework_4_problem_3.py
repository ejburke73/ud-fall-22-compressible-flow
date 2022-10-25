# Compressible Flow
# AEE 553
# Homework 4 - Problem 3
# Evan Burke

import numpy as np
from matplotlib import pyplot as plt
import shocks as ns
import isentropic as isen
import oblique as os
from homework_4_problem_2 import SimpleRamjet

M = 5
T = 217
p = 20000

theta1 = 7 #treat this as positive
delta = 1 # weak shock solution

def find_beta(M=None,gamma=1.4,theta=None):
    theta = np.deg2rad(theta)
    lamb = ((M**2-1)**2 - 3*(1 + (gamma-1)/2*M**2) * (1 + (gamma+1)/2*M**2) * np.tan(theta)**2)**0.5
    chi = ((M**2-1)**3 - 9 * (1 + (gamma-1)/2 * M**2) * (1 + (gamma-1)/2 * M**2 + (gamma+1)/4*M**4)*np.tan(theta)**2)/lamb**3
    tan_beta = (M**2 - 1 + 2*lamb*np.cos((4*np.pi*delta+np.arccos(chi))/3)) / (3 * (1 + (gamma-1)/2*M**2)*np.tan(theta))
    beta = np.arctan(tan_beta)
    beta = np.rad2deg(beta)
    print(f'Shock angle = {beta}')
    return beta

beta1 = find_beta(M=M,theta=theta1)
M1n = os.get_m1_normal(M1=M,beta=beta1)
M2n = os.get_m2_normal(M1n=M1n)
M2 = os.get_m2(M2n=M2n,beta=beta1,theta=theta1)

theta2 = 7
beta2 = find_beta(M=M2,theta=theta2)
M2np = os.get_m1_normal(M1=M2,beta=beta2)
M3n = os.get_m2_normal(M1n=M2np)
M3 = os.get_m2(M2n=M3n,beta=beta2,theta=theta2)

pt = isen.get_total_pressure(M=M,p=p)
p2_p1 = ns.get_static_pressure_ratio_normal_shock(M1=M1n)
pt2_pt1 = ns.get_total_pressure_ratio_normal_shock(M1=M1n)
p3_p2 = ns.get_static_pressure_ratio_normal_shock(M1=M2np)
pt3_pt2 = ns.get_total_pressure_ratio_normal_shock(M1=M2np)
T2_T1 = ns.get_static_temperature_ratio_normal_shock(M1=M1n)
T3_T2 = ns.get_static_temperature_ratio_normal_shock(M1=M2np)

p3_p1 = p2_p1*p3_p2
pt3_pt1 = pt3_pt2*pt2_pt1
T3 = T3_T2 * T2_T1 * T
print(f'T3 = {T3}')
p3 = p3_p1 * p
pt3 = pt3_pt2*pt2_pt1*pt
ramjet_m5 = SimpleRamjet(M1=M,theta=21,q=500,cp=1000,gamma=1.4,T=T,delta=delta)
print('\n\n')
print(f'Inlet Static Pressure Ratio:\nRamjet = {1/ramjet_m5.p1_p3}\nScramjet = {p3_p1}\n')
print(f'Inlet Total Pressure Ratio:\nRamjet = {ramjet_m5.pt3_pt1}\nScramjet = {pt3_pt1}\n')
print(f'Combustor Inlet Static Temperature:\nRamjet = {ramjet_m5.T3}\nScramjet = {T3}')

M3 = isen.get_mach_number(p_t_ratio=pt3/p3)
print(f'Combustor Mach = {M3}')