# Compressible Flow
# AEE 553
# Homework 4 - Problem 1
# Evan Burke

from cProfile import label
import numpy as np
from matplotlib import pyplot as plt
import shocks as ns
import isentropic as isen

# b
gamma = 1.4

def hugoniot(gamma=1.4,rho_ratio=None):
    p2_p1 = ( (gamma+1)/(gamma-1) * rho_ratio - 1) / ((gamma+1)/(gamma-1) - rho_ratio)
    return p2_p1

rho_ratio = np.linspace(1,5,endpoint=True)

p2_p1_h = [hugoniot(rho_ratio=r) for r in rho_ratio] # normal shock
p2_p1_i = [r**gamma for r in rho_ratio]

fig,ax = plt.subplots()
ax.plot(rho_ratio,p2_p1_h,label='Hugoniot - NS')
ax.plot(rho_ratio,p2_p1_i,label='Isentropic')
ax.legend()
ax.set_xlabel('Density Ratio')
ax.set_ylabel('Pressure Ratio')
ax.set_title('Compression vs. Density Ratio')
plt.savefig('../images/problem_1/hugoniot_vs_isentropic_compression.png')
plt.close()

machs = np.linspace(1,6,endpoint=True)
print(machs)
pt2_pt1 = [ns.get_total_pressure_ratio_normal_shock(M1=m) for m in machs] #get_total_pressure worked without a static pressure? need error handling
print(pt2_pt1)
fig,ax = plt.subplots()
ax.plot(machs,pt2_pt1)
ax.set_xlabel('Mach')
ax.set_ylabel('Total Pressure Ratio')
ax.set_title('Normal Shock Total Pressure Ratio vs. Mach Number')
plt.savefig('../images/problem_1/compression_efficiency_NS.png')
plt.close()

# d
cp = 1004.5
R = 287
ds_isen = [(cp * np.log(r**(gamma-1)) - R * np.log(pr)) for r,pr in zip(rho_ratio,p2_p1_i)]

m1 = [ns.get_upstream_mach_normal_shock(p2_p1=pr) for pr in p2_p1_h]
pt2_pt1 = [ns.get_total_pressure_ratio_normal_shock(M1=m) for m in m1]
ds_ns = [-R*np.log(ptr) for ptr in pt2_pt1]
print(ds_isen)
print(ds_ns)

fig,ax = plt.subplots()

ax.plot(rho_ratio,ds_isen,label='Isentropic')
ax.plot(rho_ratio,ds_ns,label='Normal Shock')
ax.set_title('Entropy Change vs. Density Ratio')
ax.set_xlabel('Density ratio')
ax.set_ylabel('Entropy Change [kJ/kg*K]')
ax.legend()

plt.savefig('../images/problem_1/entropy_change.png')

pr_crit = hugoniot(rho_ratio=2.5)
m1_crit = ns.get_upstream_mach_normal_shock(p2_p1=pr_crit)
print(f'Critical Mach: {m1_crit}')