#!/usr/bin/env python

from os import stat
import compressible_flow_methods as cfm

class CompressibleFlow:

    def __init__(self,M=None,u=None,a=None,p_static=None,rho_static=None,T_static=None,p_total=None,rho_total=None,T_total=None,p_ratio=None,rho_ratio=None,T_ratio=None,gamma=1.4,gas='Air',R=287):
        self.M = M
        self.u = u
        self.a = a
        self.p_static = p_static
        self.rho_static = rho_static
        self.T_static = T_static
        self.p_total = p_total
        self.rho_total = rho_total
        self.T_total = T_total
        self.p_ratio = p_ratio
        self.rho_ratio = rho_ratio
        self.T_ratio = T_ratio
        self.gamma = gamma
        self.gas = gas
        self.R = R

        while None in self.__dict__.values():
            print(self.__dict__.values())

            if not self.p_static:
                self.p_static = cfm.get_static_pressure(M=self.M,u=self.u,gamma=self.gamma,R=self.R,rho_static=self.rho_static,T_static=self.T_static,p_total=self.p_total,rho_total=self.rho_total,T_total=self.T_total,p_ratio=self.p_ratio,rho_ratio=self.rho_ratio,T_ratio=self.T_ratio)

            if not self.rho_static:
                self.rho_static = cfm.get_static_density(M=self.M,u=self.u,gamma=self.gamma,R=self.R,p_static=self.p_static,T_static=self.T_static,p_total=self.p_total,rho_total=self.rho_total,T_total=self.T_total,p_ratio=self.p_ratio,rho_ratio=self.rho_ratio,T_ratio=self.T_ratio)

            if not self.T_static:
                self.T_static = cfm.get_static_temperature(M=self.M,u=self.u,gamma=self.gamma,R=self.R,p_static=self.p_static,rho_static=self.rho_static,p_total=self.p_total,rho_total=self.rho_total,T_total=self.T_total,p_ratio=self.p_ratio,rho_ratio=self.rho_ratio,T_ratio=self.T_ratio)

            if not self.a:
                self.a = cfm.get_sonic_velocity(gamma=self.gamma,R=self.R,p=self.p_static,rho=self.rho_static,T=self.T_static)

            if not self.u:
                self.u = cfm.get_fluid_velocity(M=self.M,a=self.a,gamma=1.4,R=287,p=self.p_static,rho=self.rho_static,T=self.T_static)

            if not self.M:
                self.M = cfm.get_mach_number(u=self.u,a=self.a,gamma=self.gamma,R=self.R,p=self.p_static,rho=self.rho_static,T=self.T_static)

            if not p_total:
                self.p_total = cfm.get_total_pressure(M=self.M,u=self.u,gamma=self.gamma,R=self.R,p_static=self.p_static,rho_static=self.rho_static,T_static=self.T_static,rho_total=self.rho_total,T_total=self.T_total,p_ratio=self.p_ratio,rho_ratio=self.rho_ratio,T_ratio=self.T_ratio)

            if not T_total:
                self.T_total = cfm.get_total_temperature(M=self.M,u=self.u,gamma=self.gamma,R=self.R,p_static=self.p_static,rho_static=self.rho_static,T_static=self.T_static,p_total=self.p_total,rho_total=self.rho_total,p_ratio=self.p_ratio,rho_ratio=self.rho_ratio,T_ratio=self.T_ratio)

            if not rho_total:
                self.rho_total = cfm.get_total_density(M=self.M,u=self.u,gamma=self.gamma,R=self.R,p_static=self.p_static,rho_static=self.rho_static,T_static=self.T_static,p_total=self.p_total,T_total=self.T_total,p_ratio=self.p_ratio,rho_ratio=self.rho_ratio,T_ratio=self.T_ratio)

            if not p_ratio:
                self.p_ratio = cfm.get_pressure_ratio(M=self.M,u=self.u,gamma=self.gamma,R=self.R,p_static=self.p_static,rho_static=self.rho_static,T_static=self.T_static,p_total=self.p_total,rho_total=self.rho_total,T_total=self.T_total,rho_ratio=self.rho_ratio,T_ratio=self.T_ratio)
            
            if not T_ratio:
                self.T_ratio = cfm.get_temperature_ratio(M=self.M,u=self.u,gamma=self.gamma,R=self.R,p_static=self.p_static,rho_static=self.rho_static,T_static=self.T_static,p_total=self.p_total,rho_total=self.rho_total,T_total=self.T_total,p_ratio=self.p_ratio,rho_ratio=self.rho_ratio)

            if not rho_ratio:
                self.rho_ratio = cfm.get_density_ratio(M=self.M,u=self.u,gamma=self.gamma,R=self.R,p_static=self.p_static,rho_static=self.rho_static,T_static=self.T_static,p_total=self.p_total,rho_total=self.rho_total,T_total=self.T_total,p_ratio=self.p_ratio,T_ratio=self.T_ratio)

        print('\n\n\n')
        print('Fluid Properties:')
        print(f'Static Pressure:\t{self.p_static}\tStatic Temperature:\t{self.T_static}\tStatic Density:\t{self.rho_static}')
        print(f'Total Pressure:\t{self.p_total}\tTotal Temperature:\t{self.T_total}\tTotal Density:\t{self.rho_total}')
        print(f'Pressure Ratio:\t{self.p_ratio}\tTemperature Ratio:\t{self.T_ratio}\t Density Ratio:\t{self.rho_ratio}')
        print(f'Mach Number:\t{self.M}\tFluid Velocity:\t{self.u}\tSpeed of Sound:\t{self.a}')
        print('\n\n\n')

def get_sonic_velocity(gamma=1.4,R=287,T=None,p=None,rho=None):

    a = None
    
    while not a:

        if T:
            a = (gamma*R*T)**.5

        elif p and rho:
            a = (gamma*p/rho)**0.5
                
    print('Speed of sound = {a}')
    return a

def get_fluid_velocity(M=None,a=None,gamma=1.4,R=287,p=None,rho=None,T=None):

    u = None

    while not u:

        if M and a:
            u = M*a

        elif M and not a and T:
            u = M*(gamma*R*T)**0.5

        elif M and not T and p and rho:
            u = M*(gamma*p/rho)**0.5

    print(f'Fluid velocity = {u}')
    return u

def get_mach_number(u=None,a=None,gamma=1.4,R=287,p=None,rho=None,T=None):

    M = None

    while not M:

        if u and a:
            M = u/a

        elif u and not a and T:
            M = u/(gamma*R*T)**0.5

        elif u and not a and not T and p and rho:
            M = u/(gamma*p/rho)**0.5

    print(f'Mach number= {M}')
    return M

def get_static_pressure(M=None,u=None,gamma=1.4,R=287,rho_static=None,T_static=None,p_total=None,rho_total=None,T_total=None,p_ratio=None,rho_ratio=None,T_ratio=None):

    p_static = None

    while not p_static:

        if rho_static and T_static:
            p_static = rho_static*R*T_static
        
        elif p_total and M:
            p_static = p_total/(1 + (gamma-1) / 2 * M**2)**(gamma/(gamma-1))

        elif p_ratio and p_total:
            p_static = p_total/p_ratio

        elif p_total and T_ratio:
            p_static = p_total/(T_ratio**(gamma/(gamma-1)))

        elif p_total and T_static and T_total:
            p_static = p_total/((T_total/T_static)**(gamma/(gamma-1)))

        elif p_total and rho_ratio:
            p_static = p_total/rho_ratio**gamma

        elif p_total and rho_static and rho_total:
            p_static = p_total/(rho_total/rho_static)**gamma

        elif p_total and (u == 0):
            p_static = p_total

    print(f'Static pressure is {p_static}')
    return p_static

def get_total_pressure(M=None,u=None,gamma=1.4,R=287,p_static=None,rho_static=None,T_static=None,rho_total=None,T_total=None,p_ratio=None,rho_ratio=None,T_ratio=None):

    p_total = None

    while not p_total:

        if rho_total and T_total:
            p_total = rho_total*R*T_total

        elif p_static and M:
            p_total = (1 + (gamma-1) / 2 * M**2)**(gamma/(gamma-1))*p_static

        elif p_ratio and p_static:
            p_total = p_ratio*p_static

        elif p_static and T_ratio:
            p_total = p_static*(T_ratio**(gamma/(gamma-1)))

        elif p_static and T_static and T_total:
            p_total = p_static*((T_total/T_static)**(gamma/(gamma-1)))

        elif p_static and rho_ratio:
            p_total = p_static*rho_ratio**gamma

        elif p_static and rho_static and rho_total:
            p_total = p_static*(rho_total/rho_static)**gamma     

        elif p_static and (u == 0):
            p_total = p_static
    
    print(f'Total pressure = {p_total}')
    return p_total

def get_pressure_ratio(M=None,u=None,gamma=1.4,R=287,p_static=None,rho_static=None,T_static=None,p_total=None,rho_total=None,T_total=None,rho_ratio=None,T_ratio=None):

    p_ratio = None

    while not p_ratio:

        if M:
            p_ratio = (1 + (gamma-1) / 2 * M**2)**(gamma/(gamma-1))

        elif p_static and p_total:
            p_ratio = p_total/p_static

        elif rho_ratio:
            p_ratio = rho_ratio**(gamma)

        elif T_ratio:
            p_ratio = T_ratio**(gamma/(gamma-1))

    print(f'Pressure ratio = {p_ratio}')
    return p_ratio

def get_static_temperature(M=None,u=None,gamma=1.4,R=287,p_static=None,rho_static=None,p_total=None,rho_total=None,T_total=None,p_ratio=None,rho_ratio=None,T_ratio=None):

    T_static = None

    while not T_static:

        if p_static and rho_static:
            T_static = p_static/(rho_static*R)

        elif T_total and M:
            T_static = T_total/((1 + (gamma-1) / 2 * M**2))

        elif T_ratio and T_total:
            T_static = T_total/T_ratio

        elif p_ratio and T_total:
            T_static = T_total/p_ratio**((gamma-1)/gamma)

        elif p_static and p_total and T_total:
            T_static = T_total/(p_total/p_static)**((gamma-1)/gamma)

        elif rho_ratio and T_total:
            T_static = T_total/rho_ratio**(1/(gamma-1))

        elif rho_static and rho_total and T_total:
            T_static = T_total/(rho_total/rho_static)**(1/(gamma-1))

        elif T_total and (u == 0):
            T_static = T_total

    print(f'Static temperature = {T_static}')
    return T_static

def get_total_temperature(M=None,u=None,gamma=1.4,R=287,p_static=None,rho_static=None,T_static=None,p_total=None,rho_total=None,p_ratio=None,rho_ratio=None,T_ratio=None):

    T_total = None

    while not T_total:

        if p_total and rho_total:
            T_total = p_total/(rho_total*R)

        elif T_static and M:
            T_total = T_static*((1 + (gamma-1) / 2 * M**2))

        elif T_static and T_total:
            T_total = T_static*T_ratio

        elif p_ratio and T_static:
            T_total = T_static*p_ratio**((gamma-1)/gamma)

        elif p_static and p_total and T_static:
            T_total = T_static*(p_total/p_static)**((gamma-1)/gamma)            

        elif rho_ratio and T_static:
            T_total = T_static*rho_ratio**(1/(gamma-1))

        elif rho_static and rho_total and T_static:
            T_total = T_static*(rho_total/rho_ratio)**(1/(gamma-1))

        elif T_static and (u == 0):
            T_total = T_static

    print(f'Total temperature = {T_total}')
    return T_total

def get_temperature_ratio(M=None,u=None,gamma=1.4,R=287,p_static=None,rho_static=None,T_static=None,p_total=None,rho_total=None,T_total=None,p_ratio=None,rho_ratio=None):

    T_ratio = None

    while not T_ratio:

        if M:
            T_ratio = (1 + (gamma-1) / 2 * M**2)

        elif T_static and T_total:
            T_ratio = T_total/T_static

        elif rho_ratio:
            T_ratio = rho_ratio**(gamma-1)

        elif p_ratio:
            T_ratio = p_ratio**((gamma-1)/gamma)

    print(f'Temperature ratio = {T_ratio}')
    return T_ratio

def get_static_density(M=None,u=None,gamma=1.4,R=287,p_static=None,T_static=None,p_total=None,rho_total=None,T_total=None,p_ratio=None,rho_ratio=None,T_ratio=None):

    rho_static = None
    
    while not rho_static:

        if p_static and T_static:
            rho_static = p_static/(R*T_static)

        elif rho_total and M:
            rho_static = rho_total/((1 + (gamma-1) / 2 * M**2))**(1/(gamma-1))

        elif rho_total and rho_ratio:
            rho_static = rho_total/rho_ratio
        
        elif p_ratio and rho_total:
            rho_static = rho_total/p_ratio**(1/gamma)

        elif p_static and p_total and rho_total:
            rho_static = rho_total/(p_total/p_static)**(1/gamma)

        elif T_ratio and rho_total:
            rho_static = rho_total/T_ratio**(1/(gamma-1))

        elif T_static and T_total and rho_total:
            rho_static = rho_total/(T_total/T_static)**(1/(gamma-1))

        elif rho_total and (u == 0):
            rho_static = rho_total

    print(f'Static density = {rho_static}')
    return rho_static

def get_total_density(M=None,u=None,gamma=1.4,R=287,p_static=None,rho_static=None,T_static=None,p_total=None,T_total=None,p_ratio=None,rho_ratio=None,T_ratio=None):
    rho_total = None
    
    while not rho_total:

        if p_total and T_total:
            rho_total = p_total/(R*T_total)

        elif rho_static and M:
            rho_total = rho_static*((1 + (gamma-1) / 2 * M**2))**(1/(gamma-1))

        elif rho_static and rho_ratio:
            rho_total = rho_static*rho_ratio
        
        elif p_ratio and rho_static:
            rho_total = rho_static*p_ratio**(1/gamma)

        elif p_static and p_total and rho_static:
            rho_total = rho_static*(p_total/p_static)**(1/gamma)
           
        elif T_ratio and rho_static:
            rho_total = rho_static*T_ratio**(1/(gamma-1))

        elif T_static and T_total and rho_static:
            rho_total = rho_static*(T_total/T_static)**(1/(gamma-1))

        elif rho_static and (u == 0):
            rho_total = rho_static

    print(f'Total density = {rho_total}')
    return rho_total

def get_density_ratio(M=None,u=None,gamma=1.4,R=287,p_static=None,rho_static=None,T_static=None,p_total=None,rho_total=None,T_total=None,p_ratio=None,T_ratio=None):

    rho_ratio = None

    while not rho_ratio:

        if M:
            rho_ratio = (1 + (gamma-1) / 2 * M**2)**(1/(gamma-1))

        elif rho_static and rho_total:
            rho_ratio = rho_total/rho_static

        elif T_ratio:
            rho_ratio = T_ratio**(1/(gamma-1))

        elif p_ratio:
            rho_ratio = p_ratio**(1/gamma)

    print(f'Density ratio = {rho_ratio}')
    return rho_ratio

if __name__=='__main__':
    fluid = CompressibleFlow(M=2,p_static=101320,T_static=273)
    print(fluid.__dict__.values())
    print('I am done')