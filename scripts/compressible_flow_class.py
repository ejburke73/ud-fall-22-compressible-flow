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

if __name__=='__main__':
    fluid = CompressibleFlow(M=2,p_static=101320,T_static=273)
    print(fluid.__dict__.values())
    print('I am done')