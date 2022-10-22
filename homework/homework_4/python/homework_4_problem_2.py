# Compressible Flow
# AEE 553
# Homework 4 - Problem 2
# Evan Burke

import numpy as np
from matplotlib import pyplot as plt
import shocks as ns
import oblique as os
import isentropic as isen
from scipy.optimize import fsolve

M = 3
T = 217 # K
p = 20000 # Pa
gamma = 1.4
R = 287 # J/kg K
cp = 1000 # J/kg K
q = 500 # kJ/kg

thetas = np.linspace(0.1,34,endpoint=True)
print(thetas)

delta = 1

def find_beta(M=None,gamma=1.4,theta=None):
    theta = np.deg2rad(theta)
    lamb = ((M**2-1)**2 - 3*(1 + (gamma-1)/2*M**2) * (1 + (gamma+1)/2*M**2) * np.tan(theta)**2)**0.5
    chi = ((M**2-1)**3 - 9 * (1 + (gamma-1)/2 * M**2) * (1 + (gamma-1)/2 * M**2 + (gamma+1)/4*M**4)*np.tan(theta)**2)/lamb**3
    tan_beta = (M**2 - 1 + 2*lamb*np.cos((4*np.pi*delta+np.arccos(chi))/3)) / (3 * (1 + (gamma-1)/2*M**2)*np.tan(theta))
    beta = np.arctan(tan_beta)
    beta = np.rad2deg(beta)
    print(f'Shock angle = {beta}')
    return beta

betas = [find_beta(M=M,theta=th) for th in thetas if not np.isnan(find_beta(M=M,theta=th))]
#thetas = thetas[:-2]
print(thetas)

for beta,theta in zip(betas,thetas): # loop to solve for ideal half angle based on efficiency
    print(f'Theta = {theta}, Beta = {beta}')
    M1n = os.get_m1_normal(M1=M,beta=beta)
    M2n = os.get_m2_normal(M1n=M1n)
    M2 = os.get_m2(M2n=M2n,beta=beta,theta=theta)
    M3 = ns.get_mach_normal_shock(M1=M2)
    print('\n\n')

    
    


'''def wave_angle(gamma,M1,theta):
    func = lambda beta : np.tan(theta * np.pi/180) - 2/np.tan(beta) * (( M1 * np.sin(beta *  np.pi/180) )**2 - 1) / (M1**2 * (gamma + np.cos(2 * beta * np.pi/180)) + 2 )
    return fsolve(func,theta)

betas = [wave_angle(gamma=gamma,M1=M,theta=th*1.1) for th in thetas]

fig,ax = plt.subplots()
plt.plot(thetas,betas,'.')
plt.show()'''