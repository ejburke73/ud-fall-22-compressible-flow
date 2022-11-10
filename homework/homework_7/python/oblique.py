#!/usr/bin/env python

from re import M
import numpy as np
import shocks as ns
from scipy.optimize import fsolve

# Fully calculates post normal shock state given Mach and gamma
# Ratios can be calculated without conditions
# Post shock conditions require pre shock conditions

# Does not check for errors
# Does not solve for pre-shock mach number numerically, requires pressure ratio
# To Do: 
# should be able to only give one of the conditions and gamma and solve for everything else
# need to rename functions to be less dumb
# recast shock relations into general oblique shock relations
'''
class NormalShock:

    """
    Normal shock class needs Mach, gamma input to calculate property ratios
    Given pre-shock conditions, can get post-shock conditions
    """

    def __init__(self,M1=None,gamma=1.4,p1=None,rho1=None,T1=None,p1_t=None,T1_t=None,rho1_t=None):
        self.gamma = gamma
        self.M1 = M1
        self.M2 = get_mach_normal_shock(M1=self.M1,gamma=self.gamma)
        self.p1 = p1
        self.p2 = get_static_pressure_normal_shock(M1=self.M1,gamma=self.gamma,p1=self.p1)
        self.p2_p1 = get_static_pressure_ratio_normal_shock(M1=self.M1,gamma=self.gamma)
        self.rho1 = rho1
        self.rho2 = get_static_density_normal_shock(M1=self.M1,gamma=self.gamma,rho1=self.rho1)
        self.rho2_rho1 = get_static_density_ratio_normal_shock(M1=self.M1,gamma=self.gamma)
        self.T1 = T1
        self.T2 = get_static_temperature_normal_shock(M1=self.M1,gamma=self.gamma,T1=self.T1)
        self.T2_T1 = get_static_temperature_ratio_normal_shock(M1=self.M1,gamma=self.gamma)
        self.p1_t = p1_t
        self.p2_t = get_total_pressure_normal_shock(M1=self.M1,gamma=self.gamma,p1_t=self.p1_t)
        self.p2_t_p1_t = get_total_pressure_ratio_normal_shock(M1=self.M1,gamma=self.gamma)
        self.T1_t = T1_t
        self.T2_t = get_total_temperature_normal_shock(T1_t=self.T1_t)
        self.T2_t_T1_t = get_total_temperature_ratio_normal_shock()
        self.rho1_t = rho1_t
        self.rho2_t = get_total_density_normal_shock(p2_t_p1_t=self.p2_t_p1_t)
        self.rho2_t_rho1_t = get_total_density_ratio_normal_shock(p2_t_p1_t=self.p2_t_p1_t)
        self.sonic_area_ratio = get_sonic_area_ratio_normal_shock(p2_t_p1_t=self.p2_t_p1_t)
        self.p2_t_p1 = get_p2_total_over_p1_static_normal_shock(M2=self.M2,gamma=self.gamma,p2_p1=self.p2_p1)

class ObliqueShock:

    def __init__(self):
        pass
'''

def get_m1_normal(M1=None,beta=None):
    beta_rad = np.deg2rad(beta)
    M1n = M1*np.sin(beta_rad)
    print(f'Normal component of incoming Mach = {M1n}')
    return M1n

def get_m2_normal(M1n=None,gamma=1.4):
    M2n = ((M1n**2+2/(gamma-1)) / (2*gamma/(gamma-1)*M1n**2-1))**0.5
    print(f'Normal component of downstream Mach = {M2n}')
    return M2n

def get_m2(M2n=None,beta=None,theta=None):
    beta_rad =np.deg2rad(beta)
    theta_rad = np.deg2rad(theta)
    M2 = M2n/np.sin(beta_rad-theta_rad)
    print(f'Post oblique shock Mach number = {M2}')
    return M2

def theta_beta_mach(M=None,theta=None,beta=None,gamma=1.4):
    beta_rad = np.deg2rad(beta)
    #eq = np.tan(theta) - 1/(2*np.tan(beta)) * (M**2*np.sin(beta)**2-1)/(M**2*(gamma+np.cos(2*beta))+2)
    tanth = 2/(np.tan(beta_rad)) * (M**2*np.sin(beta_rad)**2-1)/(M**2*(gamma+np.cos(2*beta_rad))+2)
    print(f'tantheta = {tanth}')
    return tanth

def wave_angle(gamma,M1,theta):
    func = lambda beta : np.tan(theta * np.pi/180) - 2/np.tan(beta) * (( M1 * np.sin(beta *  np.pi/180) )**2 - 1) / (M1**2 * (gamma + np.cos(2 * beta * np.pi/180)) + 2 )
    return fsolve(func,53)

def get_mach_shock(M1=None,gamma=1.4,beta=90,theta=0):
    beta = np.deg2rad(beta)
    theta = np.deg2rad(theta)
    M2 = (((1 + (gamma-1)/2 * M1**2*np.sin(beta)**2) / (gamma * M1**2*np.sin(beta)**2 - (gamma-1)/2))/np.sin(beta-theta)**2)**0.5
    top = (M1**2 * np.sin(beta)**2 + 2/(gamma-1))
    bottom = (2*gamma/(gamma-1)*M1**2*np.sin(beta)**2-1)
    foo = (top/bottom)/(np.sin(beta-theta)**2)
    bar = foo**0.5
    print(f'Post normal-shock Mach number = {M2}') 
    print(foo)
    return M2

def get_upstream_mach_normal_shock(M2=None,gamma=1.4,p2_p1=None):
    M1 = ((p2_p1 + (gamma-1)/(gamma+1))*(gamma+1)/(2*gamma))**0.5
    print(f'Pre shock mach number = {M1}')
    return M1

def get_static_pressure_normal_shock(M1=None,gamma=1.4,p1=None):

    if p1:

        p2 = (1 + 2*gamma/(gamma+1)*(M1**2-1))*p1
        print(f'Post shock static pressure = {p2}')
        return p2

def get_static_pressure_ratio_normal_shock(M1=None,gamma=1.4):
    p2_p1 = (1 + 2*gamma/(gamma+1)*(M1**2-1))
    print(f'Post shock static pressure ratio = {p2_p1}')
    return p2_p1

def get_static_density_normal_shock(M1=None,gamma=1.4,rho1=None):

    if rho1:

        rho2 = (gamma+1)*M1**2/(2+(gamma-1)*M1**2)*rho1
        print(f'Post shock static density = {rho2}')
        return rho2

def get_static_density_ratio_normal_shock(M1=None,gamma=1.4):
    rho2_rho1 = (gamma+1)*M1**2/(2+(gamma-1)*M1**2)
    print(f'Post shock static density ratio = {rho2_rho1}')
    return rho2_rho1

def get_static_temperature_normal_shock(M1=None,gamma=1.4,T1=None):

    if T1:

        T2 = (1 + 2*gamma/(gamma+1) * (M1**2-1)) *  ((2 + (gamma-1) * M1**2)/((gamma+1)*M1**2))*T1
        print(f'Post shock static temperature = {T2}')
        return T2

def get_static_temperature_ratio_normal_shock(M1=None,gamma=1.4):
    T2_T1 = (1 + 2*gamma/(gamma+1) * (M1**2-1)) *  ((2 + (gamma-1) * M1**2)/((gamma+1)*M1**2))
    print(f'Post shock static temperature ratio = {T2_T1}')
    return T2_T1

def get_total_pressure_normal_shock(M1=None,gamma=1.4,p1_t=None):

    if p1_t:

        p2_t = (((gamma+1)/2*M1**2)/(1 + (gamma-1)/2 * M1**2))**(gamma/(gamma-1)) * (2*gamma/(gamma+1)*M1**2-(gamma-1)/(gamma+1))**(1/(1-gamma))*p1_t
        print(f'Post shock total pressure ratio = {p2_t}')
        return p2_t    

def get_total_pressure_ratio_normal_shock(M1=None,gamma=1.4):
    p2_t_p1_t = (((gamma+1)/2*M1**2)/(1 + (gamma-1)/2 * M1**2))**(gamma/(gamma-1)) * (2*gamma/(gamma+1)*M1**2-(gamma-1)/(gamma+1))**(1/(1-gamma))
    print(f'Post shock total pressure ratio = {p2_t_p1_t}')
    return p2_t_p1_t

def get_total_temperature_normal_shock(T1_t=None):

    if T1_t:

        T2_t = T1_t
        print(f'Post shock total temperature = {T2_t}')
        return T2_t

def get_total_temperature_ratio_normal_shock():
    T2_t_T1_t = 1
    print(f'Post shock total temperature ratio = {T2_t_T1_t}')
    return T2_t_T1_t

def get_total_density_normal_shock(p2_t_p1_t=None,rho1_t=None):

    if rho1_t:

        rho2_t = p2_t_p1_t*rho1_t
        print(f'Post shock total density = {rho2_t}')
        return rho2_t

def get_total_density_ratio_normal_shock(p2_t_p1_t=None):
    rho2_t_rho1_t = p2_t_p1_t
    print(f'Post shock total density ratio = {rho2_t_rho1_t}')
    return rho2_t_rho1_t

def get_sonic_area_ratio_normal_shock(p2_t_p1_t=None):
    sonic_area_ratio = 1/p2_t_p1_t
    print(f'Post shock sonic area ratio = {sonic_area_ratio}')
    return sonic_area_ratio

def get_p2_total_over_p1_static_normal_shock(M2=None,gamma=1.4,p2_p1=None):
    p2_t_p1 = (1 + (gamma-1)/2*M2**2)**(gamma/(gamma-1))*p2_p1
    print(f'P2t/P1 = {p2_t_p1}')
    return p2_t_p1

    '''
    # Oblique Shock Relations

    #def WaveAngle(gamma,M1,theta):
    #    func = lambda beta : np.tan(theta * np.pi/180) - 2/np.tan(beta) * (( M1 * np.sin(beta *  np.pi/180) )**2 - 1) / (M1**2 * (gamma + np.cos(2 * beta * np.pi/180)) + 2 )

    #    return fsolve(func,53)

    def M2Oblique(gamma,M1,beta,theta):
        return ( ( ( 1 + (gamma - 1) / 2 * ( M1 * np.sin(beta *  np.pi/180) )**2 ) / ( gamma * ( M1 * np.sin(beta *  np.pi/180) )**2 - (gamma - 1) / 2 ) )**0.5) / ( np.sin( (beta - theta) * np.pi/180) )

    def P2_P1Oblique(gamma,M1,beta):
        return 2 * gamma / (gamma + 1) * (M1 * np.sin(beta *  np.pi/180) )**2 - (gamma - 1) / (gamma + 1)

    def T2_T1Oblique(gamma,M1,beta):
        return (1 + (gamma - 1) / 2 * (M1 * np.sin(beta *  np.pi/180) )**2) * (2 * gamma / ( gamma - 1) * (M1 * np.sin(beta *  np.pi/180) )**2 - 1 ) / ( (M1 * np.sin(beta *  np.pi/180) )**2 * (gamma + 1)**2 / (2*(gamma-1))  )

    def Pt2_Pt1Oblique(gamma,M1,beta):
        M = M1 * np.sin(beta *  np.pi/180)
        return ( (gamma+1)/2*M**2  / ( 1 + (gamma-1)/2 * M**2) )**(gamma/(gamma-1)) * (2*gamma/(gamma+1) * M**2 - (gamma-1)/(gamma+1))**(1/(1-gamma))
    '''


if __name__=='__main__':
    #mach = get_mach_normal_shock(M1=1.5,gamma=1.4)
    #p2 = get_static_pressure_normal_shock(M1=1.5,p1=101322)
    #p2r = get_static_pressure_ratio_normal_shock(M1=1.5)
    #rho2 = get_static_density_normal_shock(M1=1.5,rho1=1.22)
    #rho2r = get_static_density_ratio_normal_shock(M1=1.5)
    #T2 = get_static_temperature_normal_shock(M1=1.5,T1=300)
    #T2r = get_static_temperature_ratio_normal_shock(M1=1.5)
    #ptr = get_total_pressure_ratio_normal_shock(M1=1.5)
    #p2t_p1 = get_p2_total_over_p1_static_normal_shock(M2=mach,p2_p1=p2r)
    #M1 = get_upstream_mach_normal_shock(M2=0.8,p2_p1=1.72)

    goo = get_mach_shock(M1=2,theta=20,beta=20.8009155)
    foo = get_mach_shock(M1=2)
    boo = get_m1_normal(M1=5,beta=29.8)
    gg = get_m2_normal(M1n=boo)
    sdf = get_m2(M2n=gg,beta=29.800,theta=20.0)
    zssdf = theta_beta_mach(M=5,beta=29.800)