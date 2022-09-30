#!/usr/bin/env python

import isentropic as isen
import shocks
import pandas as pd
from matplotlib import pyplot as plt
import numpy as np
from scipy.optimize import fsolve

# Problem 2

# Read .csv data file in as pandas dataframe, remove nan rows, convert to np array
trajectory = pd.read_csv('../HW_3_problem_2_data.csv')
trajectory = trajectory.dropna()
trajectory = trajectory.to_numpy()

# Break into specific columns
time = trajectory[:,0]
p_t_probe = trajectory[:,1]
p1 = trajectory[:,2]

# Part a
# p1 =  static pressure upstream of shock
# p_t_probe = total pressure at stagnation point = p1_t for subsonic = p2_t for supersonic

# For subsonic cases:
# p1_t/p1 = (1 + (gamma-1)/2*M1**2)**((gamma-1)/gamma)
# This is valid until the pressure ratio reaches the sonic ratio limit at M = 1

# p1_t/p* = (1 + (gamma-1)/2)**((gamma-1)/gamma) = ~1.89

# Any ratio of p_t_probe/p1 < 1.89 = subsonic
# Any ratio of p_t_probe/p1 = 1.89 = exactly sonic
# Any ratio of p_t_probe/p1 > 1.89 = subsonic --> Mach # calculated here NOT valid

# Must use:
# For supersonic cases:

# p2t/p1 =p2t/p2*p2/p1 

# p2t/p2 = (1 + (gamma-1) / 2 * M2**2)**(gamma/(gamma-1)) 
# p2_p1 = (1 + 2*gamma/(gamma+1)*(M1**2-1))
# M2**2 = ((1 + (gamma-1)/2 * M1**2) / (gamma * M1**2 - (gamma-1)/2))

# p2t/p1 = (1 + (gamma-1) / 2 * M2**2)**(gamma/(gamma-1)) * (1 + 2*gamma/(gamma+1)*(M1**2-1))
# p2t/p1 = (1 + (gamma-1) / 2 * ((1 + (gamma-1)/2 * M1**2) / (gamma * M1**2 - (gamma-1)/2)))**(gamma/(gamma-1)) * (1 + 2*gamma/(gamma+1)*(M1**2-1))
# Must be iteratively solved for M1

# p2t/p1 (M1=1) = 1.89

# Vehicle mach number must be calculated in piecewise fashion
# Generate list of pressure ratios
# All values < 1.89 -> p1_t/p1
# All values > 1.89 -> p2_t/p1
# All values < 1.89 -> p1_t/p1

ratios = [p_t/p for p_t,p in zip(p_t_probe,p1)]
sonic_ratios = isen.get_sonic_ratios()
p_t_p_star = 1/sonic_ratios[0]

def func(M1,p2_t_p1,gamma=1.4,):
    eq = (1 + (gamma-1) / 2 * ((1 + (gamma-1)/2 * M1**2) / (gamma * M1**2 - (gamma-1)/2)))**(gamma/(gamma-1)) * (1 + 2*gamma/(gamma+1)*(M1**2-1)) - p2_t_p1
    return eq

machs = [isen.get_mach_number(p_t_ratio=r) if r < p_t_p_star else float(fsolve(func,4, args=(r))) for r in ratios]
machs_simple = [isen.get_mach_number(p_t_ratio=r) for r in ratios] # If the total pressure values were valid for in front of the shock the entire time

fig,ax = plt.subplots()
plt.plot(time,machs,label='piecewise')
ax.set_xlabel('Time [s]')
ax.set_ylabel('Mach')
ax.set_title('F-16 Mach vs. Time')
plt.savefig('../images/problem_2/Mach_vs_Time_F16.png')
plt.plot(time,machs_simple,'r',label='shock-free')
ax.legend()
plt.savefig('../images/problem_2/Mach_vs_Time_F16_Piecewise_Simple.png')
plt.close()

fig,ax = plt.subplots()
plt.plot(time,ratios,label='Probe Total/Static Pressure')
plt.plot([time[0],time[-1]],[p_t_p_star,p_t_p_star],'r',label='Supersonic Line')
ax.set_xlabel('Time [s]')
ax.set_ylabel(r'$\frac{p_0}{p}$', fontsize=18)
ax.set_title('Probe Total/Static Pressure vs. Time')
ax.legend()
plt.savefig('../images/problem_2/Probe_P_ratio_vs_Time_F16.png')
plt.close()

# Part b

# Mach # experienced by probe = subsonic entire time
# Piecewise built up from aircraft Mach # when subsonic and post-normal shock Mach # when aircraft is supersonic

mach_near_probe = [M1 if M1 < 1 else shocks.get_mach_normal_shock(M1) for M1 in machs] # this is actually the mach number immediately behind the shock but not the probe because probe is stagnant

# Probe will experience same static temp as temp behind shock
# Probe will experience same total temp as temp behind shock..??? Shouldn't the static temp just be the total temp for the probe????

fig,ax = plt.subplots()
plt.plot(time,mach_near_probe,label='In Front of Probe')
ax.set_xlabel('Time [s]')
ax.set_ylabel('Mach')
ax.set_title('Region 2 Mach vs. Time')
plt.savefig('../images/problem_2/Mach_vs_Time_F16_Probe.png')
plt.plot(time,machs,'r',label='Aircraft')
ax.set_title('Region 1/2 Mach vs. Time')
ax.legend()
plt.savefig('../images/problem_2/Mach_vs_Time_F16_Probe_Aircraft.png')
plt.close()

# T = 298 at all altitudes, static temp
T_ts = [isen.get_total_temperature(T=298,M=Mi) for Mi in machs] # Total temperature in front of shock
T_ts_probe = T_ts # Total temperature does not change across a shock

T2_shock = [shocks.get_static_temperature_normal_shock(M1=Mi,T1=298) if Mi > 1 else 298 for Mi in machs] # Also the probe static temp based off of temp after normal shock 

fig,ax = plt.subplots()
plt.plot(time,T2_shock)
ax.set_xlabel('Time [s]')
ax.set_ylabel('Static Temperature [K]')
ax.set_title('Region 2 Static Temperature vs. Time')
plt.annotate('Note: Static temp is\nconstant (T=298 K)\nuntil shock forms',(0,425))
plt.savefig('../images/problem_2/T2_vs_Time_F16.png')
plt.close()

fig,ax = plt.subplots()
plt.plot(time,T_ts_probe)
ax.set_xlabel('Time [s]')
ax.set_ylabel('Probe Stagnation Temperature [K]')
ax.set_title('Probe Stagnation Temperature vs. Time')
plt.savefig('../images/problem_2/Probe_T_t_vs_Time_F16.png')
plt.close()