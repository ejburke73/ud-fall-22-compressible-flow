#!/usr/bin/env python

from cProfile import label
from webbrowser import get
import isentropic as isen
import shocks
import pandas as pd
from matplotlib import pyplot as plt
import numpy as np

# Problem 1

# Read .csv data file in as pandas dataframe, remove nan rows, convert to np array
trajectory = pd.read_csv('../HW_3_problem_2_data.csv')
trajectory = trajectory.dropna()
trajectory = trajectory.to_numpy()

# Break into specific columns
time = trajectory[:,0]
p_t = trajectory[:,1]
p = trajectory[:,2]

# Part a
# These mach numbers are UPSTREAM of the bow shock
machs = [isen.get_mach_number(p_t_ratio=p_t/p) for p_t,p in zip (p_t,p)]

fig,ax = plt.subplots()
plt.plot(time,machs)
ax.set_xlabel('Time [s]')
ax.set_ylabel('Mach')
ax.set_title('Mach vs. Time')
plt.savefig('../images/Mach_vs_Time_F16.png')
plt.close()

# Part b
mach_probe = [shocks.get_mach_normal_shock(M1=M1) for M1 in machs]


'''T_ts = [isen.get_total_temperature(M=Mi,T=298) for Mi in machs]

fig,ax = plt.subplots()
plt.plot(time,T_ts)
ax.set_xlabel('Time [s]')
ax.set_ylabel('Stagnation Temperature [K]')
ax.set_title('Stagnation Temperature vs. Time')
plt.savefig('../images/T_t_vs_Time_F16.png')
plt.close()

# This should not melt aluminum'''