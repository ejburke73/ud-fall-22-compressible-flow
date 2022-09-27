#!/usr/bin/env python

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

            if not self.p_static:
                self.p_static = self.get_static_pressure()

            if not self.rho_static:
                self.rho_static = self.get_static_density()

            if not self.T_static:
                self.T_static = self.get_static_temperature()

            if not self.a:
                self.a = self.get_sonic_velocity()

            if not self.u:
                self.u = self.get_fluid_velocity()

            if not self.M:
                self.M = self.get_mach_number()

            if not p_total:
                self.p_total = self.get_total_pressure()

            if not T_total:
                self.T_total = self.get_total_temperature()

            if not rho_total:
                self.rho_total = self.get_total_density()

            if not p_ratio:
                self.p_ratio = self.get_pressure_ratio()
            
            if not T_ratio:
                self.T_ratio = self.get_temperature_ratio()

            if not rho_ratio:
                self.rho_ratio = self.get_density_ratio()

    def get_sonic_velocity(self):

        self.a = None
        
        while not self.a:

            if self.T_static:
                self.a = (self.gamma*self.R*self.T_static)**.5

            elif self.p_static and self.rho_static:
                self.a = (self.gamma*self.p_static/self.rho_static)**0.5

            elif self.u == 0 and T_total:
                self.a = (self.gamma*self.R*self.T_total)**.5
                    
            elif self.u == 0 and self.p_total and self.rho_total:
                self.a = (self.gamma*self.p_total/self.rho_total)**0.5

        print('Speed of sound = {self.a}')
        return self.a

    def get_fluid_velocity(self):

        self.u = None

        while not self.u:

            if self.M and self.a:
                self.u = self.M*self.a

            elif self.M and not self.a and self.T_static:
                self.u = self.M*(self.gamma*self.R*self.T_static)**0.5

            elif self.M and not self.T_static and self.p_static and self.rho_static:
                self.u = self.M*(self.gamma*self.p_static/self.rho_static)**0.5

        print(f'Fluid velocity = {self.u}')
        return self.u

    def get_mach_number(self):

        self.M = None

        while not self.M:

            if self.u and self.a:
                self.M = self.u/self.a

            elif self.u and not self.a and self.T:
                self.M = self.u/(self.gamma*self.R*self.T)**0.5

            elif self.u and not self.a and not self.T and self.p and self.rho:
                self.M = self.u/(self.gamma*self.p/self.rho)**0.5

        print(f'Mach number= {self.M}')
        return self.M

    def get_static_pressure(self):

        self.p_static = None

        while not self.p_static:

            if self.rho_static and self.T_static:
                self.p_static = self.rho_static*self.R*self.T_static
            
            elif self.p_total and self.M:
                self.p_static = self.p_total/(1 + (self.gamma-1) / 2 * self.M**2)**(self.gamma/(self.gamma-1))

            elif self.p_ratio and self.p_total:
                self.p_static = self.p_total/self.p_ratio

            elif self.p_total and self.T_ratio:
                self.p_static = self.p_total/(self.T_ratio**(self.gamma/(self.gamma-1)))

            elif self.p_total and self.T_static and self.T_total:
                self.p_static = self.p_total/((self.T_total/self.T_static)**(self.gamma/(self.gamma-1)))

            elif self.p_total and self.rho_ratio:
                self.p_static = self.p_total/self.rho_ratio**self.gamma

            elif self.p_total and self.rho_static and self.rho_total:
                self.p_static = self.p_total/(self.rho_total/self.rho_static)**self.gamma

            elif self.p_total and (self.u == 0):
                self.p_static = self.p_total

        print(f'Static pressure is {self.p_static}')
        return self.p_static

    def get_total_pressure(self):

        self.p_total = None

        while not self.p_total:

            if self.rho_total and self.T_total:
                self.p_total = self.rho_total*self.R*self.T_total

            elif self.p_static and self.M:
                self.p_total = (1 + (self.gamma-1) / 2 * self.M**2)**(self.gamma/(self.gamma-1))*self.p_static

            elif self.p_ratio and self.p_static:
                self.p_total = self.p_ratio*self.p_static

            elif self.p_static and self.T_ratio:
                self.p_total = self.p_static*(self.T_ratio**(self.gamma/(self.gamma-1)))

            elif self.p_static and self.T_static and self.T_total:
                self.p_total = self.p_static*((self.T_total/self.T_static)**(self.gamma/(self.gamma-1)))

            elif self.p_static and self.rho_ratio:
                self.p_total = self.p_static*self.rho_ratio**self.gamma

            elif self.p_static and self.rho_static and self.rho_total:
                self.p_total = self.p_static*(self.rho_total/self.rho_static)**self.gamma     

            elif self.p_static and (self.u == 0):
                self.p_total = self.p_static
        
        print(f'Total pressure = {self.p_total}')
        return self.p_total

    def get_pressure_ratio(self):

        self.p_ratio = None

        while not self.p_ratio:

            if self.M:
                self.p_ratio = (1 + (self.gamma-1) / 2 * self.M**2)**(self.gamma/(self.gamma-1))

            elif self.p_static and self.p_total:
                self.p_ratio = self.p_total/self.p_static

            elif self.rho_ratio:
                self.p_ratio = self.rho_ratio**(self.gamma)

            elif self.T_ratio:
                self.p_ratio = self.T_ratio**(self.gamma/(self.gamma-1))

        print(f'Pressure ratio = {self.p_ratio}')
        return self.p_ratio

    def get_static_temperature(self):

        self.T_static = None

        while not self.T_static:

            if self.p_static and self.rho_static:
                self.T_static = self.p_static/(self.rho_static*self.R)

            elif self.T_total and self.M:
                self.T_static = self.T_total/((1 + (self.gamma-1) / 2 * self.M**2))

            elif self.T_ratio and self.T_total:
                self.T_static = self.T_total/self.T_ratio

            elif self.p_ratio and self.T_total:
                self.T_static = self.T_total/self.p_ratio**((self.gamma-1)/self.gamma)

            elif self.p_static and self.p_total and self.T_total:
                self.T_static = self.T_total/(self.p_total/self.p_static)**((self.gamma-1)/self.gamma)

            elif self.rho_ratio and self.T_total:
                self.T_static = self.T_total/self.rho_ratio**(1/(self.gamma-1))

            elif self.rho_static and self.rho_total and self.T_total:
                self.T_static = self.T_total/(self.rho_total/self.rho_static)**(1/(self.gamma-1))

            elif self.T_total and (self.u == 0):
                self.T_static = self.T_total

        print(f'Static temperature = {self.T_static}')
        return self.T_static

    def get_total_temperature(self):

        self.T_total = None

        while not self.T_total:

            if self.p_total and self.rho_total:
                self.T_total = self.p_total/(self.rho_total*self.R)

            elif self.T_static and self.M:
                self.T_total = self.T_static*((1 + (self.gamma-1) / 2 * self.M**2))

            elif self.T_static and self.T_total:
                self.T_total = self.T_static*self.T_ratio

            elif self.p_ratio and self.T_static:
                self.T_total = self.T_static*self.p_ratio**((self.gamma-1)/self.gamma)

            elif self.p_static and self.p_total and self.T_static:
                self.T_total = self.T_static*(self.p_total/self.p_static)**((self.gamma-1)/self.gamma)            

            elif self.rho_ratio and self.T_static:
                self.T_total = self.T_static*self.rho_ratio**(1/(self.gamma-1))

            elif self.rho_static and self.rho_total and self.T_static:
                self.T_total = self.T_static*(self.rho_total/self.rho_ratio)**(1/(self.gamma-1))

            elif self.T_static and (self.u == 0):
                self.T_total = self.T_static

        print(f'Total temperature = {self.T_total}')
        return self.T_total

    def get_temperature_ratio(self):

        self.T_ratio = None

        while not self.T_ratio:

            if self.M:
                self.T_ratio = (1 + (self.gamma-1) / 2 * self.M**2)

            elif self.T_static and self.T_total:
                self.T_ratio = self.T_total/self.T_static

            elif self.rho_ratio:
                self.T_ratio = self.rho_ratio**(self.gamma-1)

            elif self.p_ratio:
                self.T_ratio = self.p_ratio**((self.gamma-1)/self.gamma)

        print(f'Temperature ratio = {self.T_ratio}')
        return self.T_ratio

    def get_static_density(self):

        self.rho_static = None
        
        while not self.rho_static:

            if self.p_static and self.T_static:
                self.rho_static = self.p_static/(self.R*self.T_static)

            elif self.rho_total and self.M:
                self.rho_static = self.rho_total/((1 + (self.gamma-1) / 2 * self.M**2))**(1/(self.gamma-1))

            elif self.rho_total and self.rho_ratio:
                self.rho_static = self.rho_total/self.rho_ratio
            
            elif self.p_ratio and self.rho_total:
                self.rho_static = self.rho_total/self.p_ratio**(1/self.gamma)

            elif self.p_static and self.p_total and self.rho_total:
                self.rho_static = self.rho_total/(self.p_total/self.p_static)**(1/self.gamma)

            elif self.T_ratio and self.rho_total:
                self.rho_static = self.rho_total/self.T_ratio**(1/(self.gamma-1))

            elif self.T_static and self.T_total and self.rho_total:
                self.rho_static = self.rho_total/(self.T_total/self.T_static)**(1/(self.gamma-1))

            elif self.rho_total and (self.u == 0):
                self.rho_static = self.rho_total

        print(f'Static density = {self.rho_static}')
        return self.rho_static

    def get_total_density(self):

        self.rho_total = None
        
        while not self.rho_total:

            if self.p_total and self.T_total:
                self.rho_total = self.p_total/(self.R*self.T_total)

            elif self.rho_static and self.M:
                self.rho_total = self.rho_static*((1 + (self.gamma-1) / 2 * self.M**2))**(1/(self.gamma-1))

            elif self.rho_static and self.rho_ratio:
                self.rho_total = self.rho_static*self.rho_ratio
            
            elif self.p_ratio and self.rho_static:
                self.rho_total = self.rho_static*self.p_ratio**(1/self.gamma)

            elif self.p_static and self.p_total and self.rho_static:
                self.rho_total = self.rho_static*(self.p_total/self.p_static)**(1/self.gamma)
            
            elif self.T_ratio and self.rho_static:
                self.rho_total = self.rho_static*self.T_ratio**(1/(self.gamma-1))

            elif self.T_static and self.T_total and self.rho_static:
                self.rho_total = self.rho_static*(self.T_total/self.T_static)**(1/(self.gamma-1))

            elif self.rho_static and (self.u == 0):
                self.rho_total = self.rho_static

        print(f'Total density = {self.rho_total}')
        return self.rho_total

    def get_density_ratio(self):

        self.rho_ratio = None

        while not self.rho_ratio:

            if self.M:
                self.rho_ratio = (1 + (self.gamma-1) / 2 * self.M**2)**(1/(self.gamma-1))

            elif self.rho_static and self.rho_total:
                self.rho_ratio = self.rho_total/self.rho_static

            elif self.T_ratio:
                self.rho_ratio = self.T_ratio**(1/(self.gamma-1))

            elif self.p_ratio:
                self.rho_ratio = self.p_ratio**(1/self.gamma)

        print(f'Density ratio = {self.rho_ratio}')
        return self.rho_ratio

if __name__=='__main__':
    fluid = CompressibleFlow(M=2,p_static=101320,T_static=273)
    print(fluid.__dict__.values())
    print('I am done')