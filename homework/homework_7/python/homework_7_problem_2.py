# Compressible Flow
# AEE 553
# Homework 7 - Problem 2
# Evan Burke

import shocks as ns
import oblique as os
import isentropic as isen
import numpy as np

from scipy.integrate import odeint

gamma = 1.4

def conical_solver(M_inf=None,theta_s=None,gamma=gamma):
    # Assume a shock angle theta_s
    # Assume a freestream Mach, M_inf
    # from these, calculate flow deflection angle, delta
    delta = find_delta(M=M_inf,theta_s=theta_s,gamma=gamma)
    # using oblique shock relations (delta=theta,theta_s=beta), calculate M2
    M1n = os.get_m1_normal(M1=M_inf,beta=theta_s)
    M2n = os.get_m2_normal(M1n=M1n)
    M2 = os.get_m2(M2n=M2n,beta=theta_s,theta=delta)
    print(f'M2={M2}')
    print(f'Delta = {delta}')

    # Calculate Vr and V_th from M2 and delta
    theta_s = np.deg2rad(theta_s)
    delta = np.deg2rad(delta)

    V_p = (2/((gamma-1)*M2**2)+1)**(-0.5)
    Vp_r = V_p * np.cos(theta_s-delta)
    Vp_th = -V_p * np.sin(theta_s-delta)
    dv = [Vp_r,Vp_th]

    thetas = np.linspace(theta_s,np.deg2rad(40),num=500)

    sol = odeint(tm_ode,[Vp_r,Vp_th],thetas)
    print(thetas)

    return sol


def tm_ode(y,theta):
    v_r, v_th = y
    d2v = (v_th**2 * v_r - (gamma-1)/2 * (1 - v_r**2 - v_th**2) * (2*v_r + v_th/np.tan(theta))) / ((gamma-1)/2 * (1-v_r**2-v_th**2)-v_th**2)
    array = [v_th, d2v]
    return array

def find_delta(M=None,theta_s=None,gamma=1.4):
    theta_s = np.deg2rad(theta_s)
    tan_del = 2 / np.tan(theta_s) * (M**2 * np.sin(theta_s)**2 - 1) / (M**2 * (gamma + np.cos(2*theta_s)) + 2)
    delta = np.arctan(tan_del)
    delta = np.rad2deg(delta)
    return delta 


foo = conical_solver(M_inf=2.5,theta_s=30,gamma=gamma)
print(foo)

np.savetxt('foo.csv', foo, delimiter=',')