from signal import SIG_DFL


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
        a = sonic_velocity(gamma,R,T)
    
    M = u/a
    print(f'Mach number = {M}')
    return M

def T_total(T_static=None,M=None,gamma=1.4):
    T_total = (1 + (gamma-1) / 2 * M**2)*T_static
    print(f'Total temperature = {T_total}')
    return T_total

def p_total(p_static=None,M=None,gamma=1.4):
    p_total = (1 + (gamma-1) / 2 * M**2)**(gamma/(gamma-1))*p_static
    print(f'Total pressure = {p_total}')
    return p_total

def rho_total(rho_static=None,M=None,gamma=1.4):
    rho_total = (1 + (gamma-1) / 2 * M**2)**(1/(gamma-1))*rho_static
    print(f'Total density = {rho_total}')
    return rho_total

def T_ratio(T_static,T_total,M=None,gamma=1.4):
    T_ratio = (1 + (gamma-1) / 2 * M**2)
    print(f'T0/T = {T_ratio}')
    return T_ratio

def p_ratio(p_static=None,p_total=None,M=None,gamma=1.4):
    p_ratio = (1 + (gamma-1) / 2 * M**2)**(gamma/(gamma-1))
    print(f'p0/p = {p_ratio}')       
    return p_ratio

def rho_ratio(rho_static=None,rho_total=None,M=None,gamma=1.4):
    rho_ratio = (1 + (gamma-1) / 2 * M**2)**(1/(gamma-1))
    print(f'rho0/rho = {rho_ratio}')
    return rho_ratio

def T_static(T_total=None,M=None,gamma=1.4):
    T_static = T_total/((1 + (gamma-1) / 2 * M**2))
    print(f'Static temperature = {T_static}')
    return T_static

def p_static(p_total=None,M=None,gamma=1.4):
    p_static = p_total/(1 + (gamma-1) / 2 * M**2)**(gamma/(gamma-1))
    print(f'Static pressure = {p_static}')
    return p_static

def rho_static(rho_total=None,M=None,gamma=1.4):
    rho_static = rho_total/(1 + (gamma-1) / 2 * M**2)**(1/(gamma-1))
    print(f'Static density = {rho_static}')
    return rho_static

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

'''
def build_fluid(self):

    """
    Initialize fluid static properties
    """
    if p_static and T_static:
        rho_static = get_static_density(self)

    elif rho_static and T_static:
        p_static = get_static_pressure(self)

    elif rho_static and p_static:
        T_static = get_static_temperature(self)       

    """
    Calculate sonic velocity if not given
    """
    if not a:
        a = get_sonic_velocity(self)

    """
    Calculate Mach number if not given
    """
    if not M:
        M = get_mach_number(self)

    """
    Calculate total properties if not given
    """
    if not p_total:
        p_total = get_static_pressure(self)

    if not T_total:
        T_total = get_total_temperature(self)

    if not rho_total:
        rho_total = get_total_density(self)

    """
    Calculate ratios if not given
    """
    if not p_ratio:
        p_ratio = get_pressure_ratio(self)
    
    if not T_ratio:
        T_ratio = get_temperature_ratio(self)

    if not rho_ratio:
        rho_ratio = get_density_ratio(self)

'''

if __name__ == '__main__':
    mach = get_mach_number(u=100,a=300)
    mach = get_mach_number(u=100,T=400)
    mach = get_mach_number(u=100,rho=1.089,p=101320)
