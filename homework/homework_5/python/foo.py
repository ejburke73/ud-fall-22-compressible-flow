from cProfile import label
from cmath import pi
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.ticker import MaxNLocator
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
import shocks as ns
import oblique as os
import isentropic as isen
from scipy.optimize import fsolve

machs = np.linspace(2,10,num=9,endpoint=True)
print(machs)

gamma = 1.4

def find_theta(M=None,beta=None,gamma=1.4):
    beta = np.deg2rad(beta)
    tanth = 2 / np.tan(beta) * (M**2 * np.sin(beta)**2 - 1) / (M**2 * (gamma + np.cos(2*beta)) + 2)
    theta = np.arctan(tanth)
    theta = np.rad2deg(theta)
    #print(theta)
    return theta

data_dict = {}
data = []

for M in machs:
    prs = []
    beta_min = np.arcsin(1/M)*180/pi
    betas = np.linspace(beta_min,60,num=20,endpoint=True)

    for beta in betas:
        M1n = os.get_m1_normal(M1=M,beta=beta)
        M2n = os.get_m2_normal(M1n=M1n)
        pr = ns.get_static_pressure_ratio_normal_shock(M1=M1n)
        prs.append(pr)

    data_dict[M] = ((betas,prs))


fig,ax = plt.subplots()

for M in machs:

    data = data_dict[M]
    plt.plot(data[0],data[1],'*',label=f'Mach {M}')

ax.legend()
ax.set_xlabel('Wave Angle')
ax.set_ylabel('Pressure Ratio')
ax.set_title('Pressure Ratio vs. Wave Angle for Machs 2-10')
plt.savefig('../images/problem_1/pr_vs_beta_2D.png')


#scatter_data = [(m,b,pr) for m in machs for (b,pr) in data_dict[m]]

scatter_data = []

for m in machs:
    foo = data_dict[m]
    bs = foo[0]
    ps = foo[1]
    print(bs,ps)

    for b,p in zip(bs,ps):
        scatter_data.append((m,b,p))

print(scatter_data)

Xs = [point[0] for point in scatter_data]
Ys = [point[1] for point in scatter_data]
Zs = [point[2] for point in scatter_data]


import plotly.graph_objects as go

'''marker_data = go.Scatter3d(
    x=Xs, 
    y=Ys, 
    z=Zs, 
    marker=go.scatter3d.Marker(size=3), 
    opacity=0.8, 
    mode='markers'
)'''

data=[go.Scatter3d(x=Xs, y=Ys, z=Zs, mode='markers', marker=go.scatter3d.Marker(showscale=True), marker_color=Zs, marker_colorscale='Viridis')]

fig = go.Figure(data)

fig.update_layout(
    title='Pressure Ratio vs. Mach and Beta', 
    autosize=False,
    width=1000, 
    height=1000,
        scene=dict(
        xaxis_title='Mach',
        yaxis_title='Beta',
        zaxis_title='Pressure Ratio',
    ),
)

fig.show()
