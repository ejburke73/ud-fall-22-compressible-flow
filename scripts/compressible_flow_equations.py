#!/usr/bin/env python

class CompressibleFlow:
    def __init__(self,mach=None,u=None,t_static=None,p_static=None,rho_static=None,t_total=None,p_total=None,rho_total=None,gas='Air',R=287,gamma=1.4):
        self.mach = mach
        self.u = u
        self.t_static = t_static
        self.p_static = p_static
        self.rho_static = rho_static
        self.t_total = t_total
        self.p_total = p_total
        self.rho_total = rho_total
        self.gas = gas
        self.R = R
        self.gamma = gamma

        if self.p_static and self.rho_static and not self.t_static:
            self.t_static = self.p_static/(self.rho_static*self.R)

        if self.p_static and self.t_static and not self.rho_static:
            self.rho_static = self.p_static/(self.R*self.t_static)

        if self.rho_static and self.t_static and not self.p_static:
            self.p_static = self.rho_static*self.R*self.t_static

    def sonic_velocity(gamma=1.4,R=287,T=None):
        if T:
            a = (gamma*R*T)**.5
            print('Speed of sound = {a}')
            return a

    def ideal_gas(p=None,rho=None,T=None,R=287):
        try:
            if not p:
                p = rho*R*T
                print(f'p = {p}')
                return p
            elif not rho:
                rho = p/(R*T)
                print(f'rho = {rho}')
                return rho
            elif not T:
                T = p/(rho*R)
                print(f'T = {T}')
                return T
            else:
                print('All conditions fully defined.')
        except:
            print('Not enough inputs defined.')
        
    def mach_number(a=None,u=None,gamma=1.4,R=287,T=None):
        if not u:
            print('No velocity given!')
            return None

        if not a:
            a = CompressibleFlow.sonic_velocity(gamma,R,T)
        
        M = u/a
        print(f'Mach number = {M}')
        return M

    def T_total(T_static=None,M=None,gamma=1.4):
        T_total = (1 + (gamma-1) / 2 * M**2)*T_static
        print(f'Total temperature = {T_total}')

    def p_total(p_static=None,M=None,gamma=1.4):
        p_total = (1 + (gamma-1) / 2 * M**2)**(gamma/(gamma-1))*p_static
        print(f'Total pressure = {p_total}')

    def rho_total(rho_static=None,M=None,gamma=1.4):
        rho_total = (1 + (gamma-1) / 2 * M**2)**(1/(gamma-1))*rho_static
        print(f'Total density = {rho_total}')

    def T_ratio(T_static,T_total,M=None,gamma=1.4):
        T_ratio = (1 + (gamma-1) / 2 * M**2)
        print(f'T0/T = {T_ratio}')

    def p_ratio(p_static=None,p_total=None,M=None,gamma=1.4):
        p_ratio = (1 + (gamma-1) / 2 * M**2)**(gamma/(gamma-1))
        print(f'p0/p = {p_ratio}')       

    def rho_ratio(rho_static=None,rho_total=None,M=None,gamma=1.4):
        rho_ratio = (1 + (gamma-1) / 2 * M**2)**(1/(gamma-1))
        print(f'rho0/rho = {rho_ratio}')


if __name__ == '__main__':
    gas = CompressibleFlow(p_static=101300,t_static=287)
    CompressibleFlow.ideal_gas(p=101300,rho=1.22)
    CompressibleFlow.mach_number(a=1000,u=432)
    CompressibleFlow.mach_number(T=422,u=432)