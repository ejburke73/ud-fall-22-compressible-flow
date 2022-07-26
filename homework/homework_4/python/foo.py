M = 3
T = 217 # K
p = 20000 # Pa
gamma = 1.4
R = 287 # J/kg K
cp = 1000 # J/kg K
q = 500 # kJ/kg

thetas = np.linspace(1,34,num=67,endpoint=True)
print(thetas)

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

betas = [find_beta(M=M,theta=th) for th in thetas if not np.isnan(find_beta(M=M,theta=th))]

fig,ax = plt.subplots()
plt.plot(thetas,betas)
plt.show()

def rayleigh(M2,M1,Tt2_Tt1,gamma=1.4):
    eq = (M2/M1)**2 * ((1 + gamma*M1**2)/(1 + gamma*M2**2))**2 * ((1 + (gamma-1)/2*M2**2)/(1 + (gamma-1)/2 * M1**2)) - Tt2_Tt1
    return eq

efficiency = []

for beta,theta in zip(betas,thetas): # loop to solve for ideal half angle based on efficiency
    print(f'Theta = {theta}, Beta = {beta}')
    M1n = os.get_m1_normal(M1=M,beta=beta)
    M2n = os.get_m2_normal(M1n=M1n)
    p2_p1 = ns.get_static_pressure_ratio_normal_shock(M1=M1n)
    pt2_pt1 = ns.get_total_pressure_ratio_normal_shock(M1=M1n)
    M2 = os.get_m2(M2n=M2n,beta=beta,theta=theta)
    M3 = ns.get_mach_normal_shock(M1=M2)
    p3_p2 = ns.get_static_pressure_ratio_normal_shock(M1=M2)
    pt3_pt2 = ns.get_total_pressure_ratio_normal_shock(M1=M2)
    p1_p3 = 1/p2_p1 * 1/p3_p2
    pt3_pt1 = pt2_pt1 * pt3_pt2
    T2_T1 = ns.get_static_temperature_ratio_normal_shock(M1=M1n)
    T3_T2 = ns.get_static_temperature_ratio_normal_shock(M1=M2)
    T3_T1 = T3_T2 * T2_T1
    print(f'T3_T1 = {T3_T1}')
    T3 = T3_T2 * T2_T1 * T
    print(f'T3 = {T3}')
    Tt3 = isen.get_total_temperature(M=M3,T=T3)
    Tt4 = Tt3 + (1000*q)/cp
    print(f'Tt4 = {Tt4}')
    Tt4_Tt3 = Tt4/Tt3
    print(f'Tt4/Tt3 = {Tt4_Tt3}')
    # The commmented lines are for solving with Rayleigh
    #M4 = float(fsolve(rayleigh,.5,args=(M3,Tt4/Tt3)))
    #print(f'M4 = {M4}')
    T4 = float(isen.get_static_temperature(M=M3,T_t=Tt4))
    eta = 1 - ((p1_p3)**((gamma-1)/gamma) * (T4 - (pt3_pt1)**((gamma-1)/gamma)*T3) / (T4-T3))
    print(f'Scramjet efficiency = {eta}')
    efficiency.append(float(eta))
    print('\n\n')