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
mach_probe = [shocks.get_mach_normal_shock(M1=Mi) if Mi > 1 else Mi for Mi in machs]
print(mach_probe)

fig,ax = plt.subplots()
plt.plot(time,mach_probe)
ax.set_xlabel('Time [s]')
ax.set_ylabel('Mach')
ax.set_title('Probe Mach vs. Time')
plt.savefig('../images/Probe_Mach_vs_Time_F16.png')
plt.close()

# T = 298 at all altitudes, static temp
T_ts = [isen.get_total_temperature(T=298,M=Mi) for Mi in machs] # Total temperature in front of shock
#T_ts_probe = T_ts # Total temperature at probe behind shock, does not change because shock is adiabatic
T_ts_probe = [isen.get_total_temperature(M=Mi,T=298) if Mi <1 else shocks.get_total_temperature_normal_shock(T1_t=T_ts) for Mi,T_ts in zip(machs,T_ts)]
T_probe = [isen.get_static_temperature(M=Mi,T_t=T_t_p) if Mi < 1 else shocks.get_static_temperature_normal_shock(M1=Mi,T1=298) for Mi,T_t_p in zip(machs,T_ts_probe)]
# Assuming that we have stagnation conditions at prove tip -> total temp of probe = static temp of probe = total temp before shock
print(T_probe)


fig,ax = plt.subplots()
plt.plot(time,T_ts_probe)
ax.set_xlabel('Time [s]')
ax.set_ylabel('Probe Stagnation Temperature [K]')
ax.set_title('Probe Stagnation Temperature vs. Time')
plt.savefig('../images/Probe_T_t_vs_Time_F16.png')
#plt.show()
plt.close()

# I struggle with this one because T_static is not increasing while subsonic

fig,ax = plt.subplots()
plt.plot(time,T_probe)
ax.set_xlabel('Time [s]')
ax.set_ylabel('Probe Static Temperature [K]')
ax.set_title('Probe Static Temperature vs. Time')
plt.savefig('../images/Probe_T_vs_Time_F16.png')
#plt.show()

for foo,bar in zip(T_ts,T_ts_probe):
    print(foo-bar)

'''fig,ax = plt.subplots()
#plt.plot(time,T_ts_probe)
plt.plot(time,T_probe)
#plt.plot(time,T_ts,'+b')
ax.set_xlabel('Time [s]')
ax.set_ylabel('Stagnation Temperature [K]')
ax.set_title('Static Temperature vs. Time')
#ax2=ax.twinx()
#ax2.plot(time,machs,'r')
#plt.savefig('../images/T_t_vs_Time_F16.png')
plt.show()
plt.close()


fig,ax = plt.subplots()
plt.plot(time,T_ts)
ax.set_xlabel('Time [s]')
ax.set_ylabel('Stagnation Temperature [K]')
ax.set_title('Stagnation Temperature vs. Time')
plt.savefig('../images/T_t_vs_Time_F16.png')
plt.close()'''

# This should not melt aluminum'''