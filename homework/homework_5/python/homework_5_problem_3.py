# Compressible Flow
# AEE 553
# Homework 5 - Problem 3
# Evan Burke


from sympy import solve, nsolve, symbols, pprint, asin, atan, cot, sin, cos, tan
import numpy as np
from matplotlib import pyplot as plt
import shocks as ns
import oblique as os
import isentropic as isen
from scipy.optimize import fsolve

M1 = 3
p1 = 1 # atm
th2 = 20 # deg
th3 = -15 # deg

gamma = 1.4
delta = 1 # weak shock solution

def find_beta(M=None,gamma=1.4,theta=None,delta=1):
    theta = np.deg2rad(theta)
    lamb = ((M**2-1)**2 - 3*(1 + (gamma-1)/2*M**2) * (1 + (gamma+1)/2*M**2) * np.tan(theta)**2)**0.5
    chi = ((M**2-1)**3 - 9 * (1 + (gamma-1)/2 * M**2) * (1 + (gamma-1)/2 * M**2 + (gamma+1)/4*M**4)*np.tan(theta)**2)/lamb**3
    tan_beta = (M**2 - 1 + 2*lamb*np.cos((4*np.pi*delta+np.arccos(chi))/3)) / (3 * (1 + (gamma-1)/2*M**2)*np.tan(theta))
    beta = np.arctan(tan_beta)
    beta = np.rad2deg(beta)
    print(f'Shock angle = {beta}')
    return beta

beta2 = find_beta(M=M1,theta=th2)
beta3 = find_beta(M=M1,theta=abs(th3))

M1n2 = os.get_m1_normal(M1=M1,beta=beta2)
M1n3 = os.get_m1_normal(M1=M1,beta=beta3)

M2n = os.get_m2_normal(M1n=M1n2)
M3n = os.get_m2_normal(M1n=M1n3)

M2 = os.get_m2(M2n=M2n,beta=beta2,theta=th2)
M3 = os.get_m2(M2n=M3n,beta=beta3,theta=abs(th3))

p2_p1 = ns.get_static_pressure_ratio_normal_shock(M1=M1n2)
p3_p1 = ns.get_static_pressure_ratio_normal_shock(M1=M1n3)

p2 = p2_p1  * p1
p3 = p3_p1  * p1

# th3 + th4 = th2 + th4p
# p4 = p4p

th2 = np.deg2rad(th2)
th3 = np.deg2rad(th3)

th4, th4p, b4, b4p, p4 = symbols('th4,th4p,b4,b4p,p4')

eq_1 = p4/p3 - 1 + 2*gamma/(gamma+1) * ((M3*sin(b4))**2-1)
eq_2 = p4/p2 - 1 + 2*gamma/(gamma+1) * ((M2*sin(b4p))**2-1)
eq_3 = tan(th4) - 2*cot(b4) * (M3**2*sin(b4)**2-1) / (M3**2 *(gamma + cos(2*b4)) + 2); 
eq_4 = tan(th4p) - 2*cot(b4p) * (M2**2*sin(b4p)**2-1) / (M2**2 *(gamma + cos(2*b4p)) + 2); 

eq_5 = solve(eq_1,b4)
eq_6 = solve(eq_2,b4p)
eq_7 = solve(eq_3,th4)
eq_8 = solve(eq_4,th4p)
print(type(eq_1))
foo = eq_1.subs(p4,8.35)
print(foo)
print(solve(foo,b4))

'''th4, th4p = symbols('th4,th4p')

theta_eq = th4 - th4p + th3 - th2
pprint(theta_eq)

beta4, beta4p, p4 = symbols('beta4,beta4p,p4')

b4_eq = beta4 - asin( ( ((p4/p3-1)*(gamma+1)/(2*gamma)+1)/M3**2 )**0.5 )
print('\n\n')
pprint(b4_eq)
b4p_eq = beta4p - asin( ( ((p4/p2-1)*(gamma+1)/(2*gamma)+1)/M2**2 )**0.5 )
print('\n\n')
pprint(b4p_eq)

th4_eq = th4 - atan(2*cot(beta4) * ((M3**2*sin(beta4)**2-1)/(M3**2*(gamma+cos(2*beta4))+2)) )
print('\n\n')
pprint(th4_eq)

th4p_eq = th4p - atan(2*cot(beta4p) * ((M2**2*sin(beta4p)**2-1)/(M2**2*(gamma+cos(2*beta4p))+2)) )

print('\n\n')
pprint(th4p_eq)

th4,th4p,beta4,beta4p,p4 = nsolve((theta_eq,th4_eq,th4p_eq,b4_eq,b4p_eq),(th4,th4p,beta4,beta4p,p4),(5,5,.25,.25,4))

type(th4)'''
#th4 = np.rad2deg(th4)

'''print(f'Theta4 = {th4*180/np.pi}')
print(f'Theta4p = {th4p*180/np.pi}')
print(f'Beta4 = {beta4*180/np.pi}')
print(f'Beta4p = {beta4p*180/np.pi}')
print(f'p4 = {p4}')

print(f'M2={M2}')

print(f'M2 = {M2}, M3={M3}')
print(f'p2 = {p2}, p3={p3}')'''