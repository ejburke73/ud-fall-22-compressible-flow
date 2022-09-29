#!/usr/bin/env python

from cProfile import label
from webbrowser import get
import isentropic as isen
import pandas as pd
from matplotlib import pyplot as plt
import numpy as np

# Problem 1

# Read .csv data file in as pandas dataframe, remove nan rows, convert to np array
trajectory = pd.read_csv('../HW_3_problem_1_data.csv',skiprows=[1])
trajectory = trajectory.dropna()
trajectory = trajectory.to_numpy()

traj_dict = {'time' : 0, 'h' : 1, 'u' : 2, 'T' : 3, 'p' : 4, 'rho' : 5}


# Part A
mach_isen = [isen.get_mach_number(u=u,T=T,p=p) for u, T, p in zip (trajectory[:,traj_dict['u']],trajectory[:,traj_dict['T']],trajectory[:,traj_dict['p']])]
mach_roomtemp = [isen.get_mach_number(u=u,T=298,p=p) for u,p in zip(trajectory[:,traj_dict['u']],trajectory[:,traj_dict['p']])]
mach_diff = [abs(Mi-Mr) for Mi,Mr in zip(mach_isen,mach_roomtemp)]
percent_diff = [diff/mi*100 for diff,mi in zip(mach_diff,mach_isen)]

fig, ax = plt.subplots()
plt.plot(trajectory[:,traj_dict['time']],mach_isen, label = 'Isentropic Mach')
plt.plot(trajectory[:,traj_dict['time']],mach_roomtemp, label = 'Room Temp Mach')
plt.plot(trajectory[:,traj_dict['time']],mach_diff, label =  'Absolute Difference')
ax.set_xlabel('Time [s]')
ax.set_ylabel('Mach')
ax.set_title('Mach vs. Time')
ax.legend()
plt.savefig('../images/Mach_vs_Time.png')
plt.close()

fig, ax = plt.subplots()
plt.plot(trajectory[:,traj_dict['time']],percent_diff)
ax.set_xlabel('Time [s]')
ax.set_ylabel('Mach Error [%]')
ax.set_title('Mach Error vs. Time')
plt.savefig('../images/Mach_Error_vs_Time.png')
plt.close()

# Part B
# Find the Mach number where rho_inf = rho_0 -> I know it's M = 0
mach_rho_equal = isen.get_mach_number(rho_t_ratio=1) 

# Calculate time history of rho_inf/rho_t
# isentropic function returns rho_t/rho_inf, take inverse
rho_rho_t = [1/isen.get_density_ratio(M=M) for M in mach_isen]

# Defining a flow as compressible when density changes > 5% 
# Need time in flight when Mach number when rho_inf/rho_t = 0.95
mach_rho_095 = isen.get_mach_number(rho_t_ratio=1/0.95)

closest_mach_095 = min(mach_isen, key=lambda x:abs(x-mach_rho_095))
print(f'Mach at closest time = {closest_mach_095}')
idx_095 = mach_isen.index(closest_mach_095)

closest_time_095 = trajectory[idx_095,traj_dict['time']]
print(f'Closest Time = {closest_time_095}')

closest_rho_ratio_095 = rho_rho_t[idx_095]
print(f'rho_inf/rho_t at closest time = {closest_rho_ratio_095}')


fig,ax = plt.subplots()
plt.plot(trajectory[:,traj_dict['time']],rho_rho_t)
ax.set_xlabel('Time [s]')
ax.set_ylabel(r'$\frac{\rho_\infty}{\rho_0}$', fontsize=18)
ax.set_title(r'$\frac{\rho_\infty}{\rho_0}$ vs Time')
plt.savefig('../images/rho_rho_t_vs_Time.png')
plt.close()

fig,ax = plt.subplots()
plt.plot(mach_isen,rho_rho_t)
ax.set_xlabel('Mach')
ax.set_ylabel(r'$\frac{\rho_\infty}{\rho_0}$', fontsize=18)
ax.set_title(r'$\frac{\rho_\infty}{\rho_0}$ vs Mach')
plt.savefig('../images/rho_rho_t_vs_Mach.png')
plt.close()

fig,ax1 = plt.subplots()
ax2 = ax1.twinx()
ax1.plot(trajectory[:,traj_dict['time']],rho_rho_t)
ax2.plot(trajectory[:,traj_dict['time']],mach_isen,'r')
ax1.set_xlabel('Time [s]')
ax1.set_ylabel(r'$\frac{\rho_\infty}{\rho_0}$', fontsize=18)
ax2.set_ylabel('Mach')
ax1.set_title(r'$\frac{\rho_\infty}{\rho_0}$ vs Time')
plt.savefig('../images/rho_rho_t_and_Mach_vs_Time.png')
plt.close()
ax1.plot(closest_time_095,closest_rho_ratio_095,'kd')
ax2.plot(closest_time_095,closest_mach_095,'bo')
ann1 = ax1.annotate(f'rho_inf/rho_t = {np.round(closest_rho_ratio_095,2)}\nt = {np.round(closest_time_095)}',(closest_time_095+5,closest_rho_ratio_095))
ann2 = ax2.annotate(f'M = {np.round(closest_mach_095,2)}\nt = {np.round(closest_time_095)}',(closest_time_095+5,0))
plt.savefig('../images/rho_rho_t_and_Mach_vs_Time_marked.png')
plt.close()
ann1.remove()
ann2.remove()
ann1 = ax1.annotate(f'rho_inf/rho_t = {np.round(closest_rho_ratio_095,2)}\nt = {np.round(closest_time_095)}',(closest_time_095+1,closest_rho_ratio_095))
ann2 = ax2.annotate(f'M = {np.round(closest_mach_095,2)}\nt = {np.round(closest_time_095)}',(closest_time_095+1,0))
ax1.set_xlim(left=0,right=10)
ax2.set_xlim(left=0,right=10)
plt.savefig('../images/rho_rho_t_and_Mach_vs_Time_marked_zoom.png')
plt.close()

p_t_isen = [isen.get_total_pressure(M=Mi,p=trajectory[i,traj_dict['p']])/1000 for i,Mi in enumerate(mach_isen)]
p_t_bernoulli = [(p + 0.5*rho*u**2)/1000 for p,rho,u in zip(trajectory[:,traj_dict['p']],trajectory[:,traj_dict['rho']],trajectory[:,traj_dict['u']])]

fig,ax = plt.subplots()
plt.plot(trajectory[:,traj_dict['time']],p_t_isen,label='Isentropic Total Pressure')
plt.plot(trajectory[:,traj_dict['time']],p_t_bernoulli,label='Bernoulli Total Pressure')
ax.set_xlabel('Time [s]')
ax.set_ylabel(r'$p_0$ [kPa]')
ax.set_title('Total Pressure vs. Time')
ax.legend()
plt.savefig('Pt_vs_Time_Isen_Bernoulli.png')
plt.close()

#  Part C
mach_target = 6

close_machs = [Mi for Mi in mach_isen if Mi > 0.99*mach_target and Mi < 1.01*mach_target]
print(f'Number of points within +- 1% of target mach = {len(close_machs)}')
print(f'Machs within +- 1% of target mach = {close_machs}')

closest_mach_M6 = min(mach_isen, key=lambda x:abs(x-mach_target))
print(f'Closest Mach to Mach 6 = {closest_mach_M6}')
idx_M6 = mach_isen.index(closest_mach_M6)
closest_time_M6 = trajectory[idx_M6,traj_dict['time']]
print(f'Closest time to Mach 6 = {closest_time_M6}')
closest_u_M6 = trajectory[idx_M6,traj_dict['u']]
print(f'Closest freestream velocity to Mach 6 = {closest_u_M6}')
closest_T_M6 = trajectory[idx_M6,traj_dict['T']]
T_t_M6 = isen.get_total_temperature(M=closest_mach_M6,T=closest_T_M6)
T_t_T_M6 = isen.get_temperature_ratio(M=closest_mach_M6)

# This does not seem like a feasible temperature for the blankets to generate

# If this were achievable, the static temperature in the nozzle throat:
T_sonic_blanket = isen.get_static_temperature(M=1,T_t=T_t_M6)

# Part D
# Assuming matching M and u in the tunnel
# Assuming same size vehicle
# Assuming same Reynolds number
# Are surface temps the same?

# I think so