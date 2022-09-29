#!/usr/bin/env python

#from shocks import NormalShock

# Does not check for compatibility between states
# Does not catch all errors

# To Do:
# Write unit tests

class CompressibleFlow:
    r"""CompressibleFlow class
    
    This class solves for all possible thermodynamic states of a compressible
    flow using isentropic relations. These relations assume the working fluid
    is a calorically perfect gas (CPG) with constant specific heats, :math:`c_p`
    and :math:`c_v`. Assumptions also used in the derivation of these relationships:
    steady flow, inviscid, no heat transfer/adiabatic, no irreversible work
    modes, and 1-D flow (all properties only vary in x-direction).
    The flow itself does not need to be isentropic for the methods
    utilized to be valid -- the total conditions are idealized representations
    of what would happen if the flow were isentropically brought to a rest.
    If the flow is isentropic, direct use of total conditions betweeen two
    points on a streamline is valid. Use of isentropic relations in non-
    isentropic flows can still be useful for understanding the bulk trends
    of a compressible flow problem.
    """

    def __init__(self,M=None,u=None,a=None,p=None,rho=None,T=None,p_t=None,
                rho_t=None,T_t=None,p_t_ratio=None,rho_t_ratio=None,
                T_t_ratio=None,gamma=1.4,gas='Air',R=287):
        """
        Parameters
        ----------
        M : float or int, optional
            Mach number, fluid velocity/speed of sound, by default None
        u : float or int or None, optional
            fluid velocity, length/time, by default None
        a : float or int or None, optional
            fluid speed of sound, by default None
        p : float or int or None, optional
            fluid static pressure, force/area, by default None
        rho : float or int or None, optional
            fluid static density, mass/volume, by default None
        T : float or int or None, optional
            fluid static temperature, degrees, by default None
        p_t : float or int or None, optional
            fluid total pressure, force/area, by default None
        rho_t : float or int or None, optional
            fluid total density, mass/volume, by default None
        T_t : float or int or None, optional
            fluid total temperature, degrees, by default None
        p_t_ratio : float or int or None, optional
            ratio of total pressure to static pressure, by default None
        rho_t_ratio :float or int or None, optional
            ratio of total density to static density, by default None
        T_t_ratio :float or int or None, optional
            ratio of total temperature to static temperature, by default None
        gamma : float
            ratio of specific heats, cp/cv, by default 1.4 for calorically perfect air
        gas : str, optional
            working fluid, eventually to be updated with other gas library, by default 'Air'
        R : int or float, optional
            specific gas constant, by default 287 J/kg*K for calorically perfect air
        """
        self.M = M
        self.u = u
        self.a = a
        self.p = p
        self.rho = rho
        self.T = T
        self.p_t = p_t
        self.rho_t = rho_t
        self.T_t = T_t
        self.p_t_ratio = p_t_ratio
        self.rho_t_ratio = rho_t_ratio
        self.T_t_ratio = T_t_ratio
        self.gamma = gamma
        self.gas = gas
        self.R = R

    def isentropic_state(self):
        """Calculates thermodynamic state to max extent possible.
        
        Depending on arguments given when instantiating CompressibleFlow object,
        different end states are viable for defining the object. Without fully
        defining the thermodynamic state (i.e., providing two of :math:`p, 
        \\rho, T`), only the property ratios are able to be calculated, assuming 
        a Mach number is given. This function determines which method of
        calculation is appropriate for the given inputs, and performs the
        appropriate calculations.
        """
        mode = self.__isentropic_input_check()
        print(f'Calculation mode: {mode}')

        if mode == 0:
            self.__isentropic_state_method_0()
                
        elif mode == 1:
            self.__isentropic_state_method_1()

        elif mode == 2:
            self.__isentropic_state_method_2()

    def __isentropic_state_method_0(self):
        """Calculates total to static ratios for :math:`p, \\rho, T` based on
        input Mach and :math:`\\gamma`.
        """
        self.p_t_ratio = get_pressure_ratio(M=self.M,gamma=self.gamma,R=self.R)
        self.rho_t_ratio = get_density_ratio(M=self.M,gamma=self.gamma,R=self.R)
        self.T_t_ratio = get_temperature_ratio(M=self.M,gamma=self.gamma,R=self.R)
        
    def __isentropic_state_method_1(self):
        """Calculates all thermodynamic properties based on inputs.
        Requires two of :math:`p, \\rho, T` and one of :math:`M, u`.
        """
        while None in self.__dict__.values():
            if not self.p:
                self.p = get_static_pressure(gamma=self.gamma,R=self.R,M=self.M,u=self.u,rho=self.rho,T=self.T,p_t=self.p_t,rho_t=self.rho_t,T_t=self.T_t,p_t_ratio=self.p_t_ratio,rho_t_ratio=self.rho_t_ratio,T_t_ratio=self.T_t_ratio)
            if not self.rho:
                self.rho = get_static_density(gamma=self.gamma,R=self.R,M=self.M,u=self.u,p=self.p,T=self.T,p_t=self.p_t,rho_t=self.rho_t,T_t=self.T_t,p_t_ratio=self.p_t_ratio,rho_t_ratio=self.rho_t_ratio,T_t_ratio=self.T_t_ratio)
            if not self.T:
                self.T = get_static_temperature(gamma=self.gamma,R=self.R,M=self.M,u=self.u,p=self.p,rho=self.rho,p_t=self.p_t,rho_t=self.rho_t,T_t=self.T_t,p_t_ratio=self.p_t_ratio,rho_t_ratio=self.rho_t_ratio,T_t_ratio=self.T_t_ratio)
            if not self.a:
                self.a = get_sonic_velocity(gamma=1.4,R=287,u=self.u,p=self.p,rho=self.rho,T=self.T,p_t=self.p_t,rho_t=self.rho_t,T_t=self.T_t)
            if not self.M:
                self.M = get_mach_number(gamma=1.4,R=287,a=self.a,u=self.u,p=self.p,rho=self.rho,T=self.T,p_t_ratio=self.p_t_ratio,rho_t_ratio=self.rho_t_ratio,T_t_ratio=self.T_t_ratio)
            if not self.u:
                self.u = get_fluid_velocity(gamma=1.4,R=287,M=self.M,a=self.a,p=self.p,rho=self.rho,T=self.T)
            if not self.p_t:
                self.p_t = get_total_pressure(gamma=self.gamma,R=self.R,M=self.M,u=self.u,p=self.p,rho=self.rho,T=self.T,rho_t=self.rho_t,T_t=self.T_t,p_t_ratio=self.p_t_ratio,rho_t_ratio=self.rho_t_ratio,T_t_ratio=self.T_t_ratio)
            if not self.rho_t:
                self.rho_t = get_total_density(gamma=self.gamma,R=self.R,M=self.M,u=self.u,p=self.p,rho=self.rho,T=self.T,p_t=self.p_t,T_t=self.T_t,p_t_ratio=self.p_t_ratio,rho_t_ratio=self.rho_t_ratio,T_t_ratio=self.T_t_ratio)
            if not self.T_t:
                self.T_t = get_total_temperature(gamma=self.gamma,R=self.R,M=self.M,u=self.u,p=self.p,rho=self.rho,T=self.T,p_t=self.p_t,rho_t=self.rho_t,p_t_ratio=self.p_t_ratio,rho_t_ratio=self.rho_t_ratio,T_t_ratio=self.T_t_ratio)
            if not self.p_t_ratio:
                self.p_t_ratio = get_pressure_ratio(gamma=self.gamma,R=self.R,M=self.M,u=self.u,p=self.p,rho=self.rho,T=self.T,p_t=self.p_t,rho_t=self.rho_t,T_t=self.T_t,rho_t_ratio=self.rho_t_ratio,T_t_ratio=self.T_t_ratio)
            if not self.rho_t_ratio:
                self.rho_t_ratio = get_density_ratio(gamma=self.gamma,R=self.R,M=self.M,u=self.u,p=self.p,rho=self.rho,T=self.T,p_t=self.p_t,rho_t=self.rho_t,T_t=self.T_t,p_t_ratio=self.p_t_ratio,T_t_ratio=self.T_t_ratio)
            if not self.T_t_ratio:
                self.T_t_ratio = get_temperature_ratio(gamma=self.gamma,R=self.R,M=self.M,u=self.u,p=self.p,rho=self.rho,T=self.T,p_t=self.p_t,rho_t=self.rho_t,T_t=self.T_t,p_t_ratio=self.p_t_ratio,rho_t_ratio=self.rho_t_ratio)

    def __isentropic_state_method_2(self):
        """Calculates all total to static ratios and Mach number based on
        input of one total to static ratio and gamma.
        """
        while None in [self.M,self.p_t_ratio,self.rho_t_ratio,self.T_t_ratio]:
            if not self.M:
                self.M = get_mach_number(gamma=self.gamma,p_t_ratio=self.p_t_ratio,rho_t_ratio=self.rho_t_ratio,T_t_ratio=self.T_t_ratio)
            if not self.p_t_ratio:
                self.p_t_ratio = get_pressure_ratio(M=self.M,gamma=self.gamma)
            if not self.rho_t_ratio:
                self.rho_t_ratio = get_density_ratio(M=self.M,gamma=self.gamma)
            if not self.T_t_ratio:
                self.T_t_ratio = get_temperature_ratio(M=self.M,gamma=self.gamma)
            pass

    def __isentropic_input_check(self):
        """Checks inputs to determine which method of calculations are performed.
        
        This function determines which method of calculation is appropriate for 
        the given inputs, and outputs an integer variable associated with said
        method.

        Returns
        -------
        mode : int
            Integer describing which method of calculation to use to
            evaluate state
        """
        mode = None

        if self.M:

            if not any([self.p,self.rho,self.T]):
                mode = 0
                print('Ideal gas not fully defined -- calculating thermodynamic ratios only.')
                return mode
            
            elif sum(1 for i in [self.p,self.rho,self.T] if i is not None)>1:
                mode = 1
                print('Ideal gas fully defined -- calculating all values.')
                return mode

        elif sum(1 for i in [self.p_t_ratio,self.rho_t_ratio,self.T_t_ratio] if i is not None)>0:
            mode = 2
            print('No Mach number defined -- using given total ratio to calculate other thermodynamic ratios.')
            return mode

        elif not self.M:
            raise ValueError('No Mach number defined, ideal gas not fully defined, no total ratios defined -- exiting now.')
        
    def report_outputs(self):
        for key,value in zip(self.__dict__.keys(),self.__dict__.values()):
            print(str(key) + ': ' + str(value))

    def shock(self):
        #shocked = NormalShock(M1=self.M,gamma=self.gamma,p1_static=self.p,rho1_static=self.rho_static,T1_static=self.T_static,p1_total=self.p_total,T1_total=self.T_total,rho1_total=self.rho_total)
        #post_shock = CompressibleFlow(M=shocked.M2,p=shocked.p2_static,rho_static=shocked.rho2_static,T_static=shocked.T2_static)
        #return shocked, post_shock
        pass

def get_sonic_velocity(gamma=1.4,R=287,u=None,p=None,rho=None,T=None,p_t=None,rho_t=None,T_t=None):
    """Calculates the speed of sound, :math:`a`, for a given fluid.
    Defaults to air with \gamma = 1.4, R = 287 J/kg*K.

    This function can solve for the speed of sound in several ways.

    Given :math:`T`:

        * :math:`a = \\sqrt{\\gamma*R*T}`.

    Given :math:`P` and :math:`\\rho` 
        
        * :math:`a = \\sqrt{\\gamma*p/\\rho}`.

    If the fluid velocity or Mach number are 0 the flow is stagnant and
    the speed of sound in the fluid is calculated using the above relations
    with total conditions:

        * :math:`a = \\sqrt{\\gamma*R*T_t}`

        * :math:`a = \\sqrt{\\gamma*p_t/\\rho_t}`

    Parameters
    ----------
    gamma : float, optional
        *See* :paramref:`.CompressibleFlow.gamma`, by default 1.4 for 
        calorically perfect air
    R : int or float
        *See* :paramref:`.CompressibleFlow.r`, by default 287 J/kg*K
    u: float or int or None, optional
        *See* :paramref:`.CompressibleFlow.u`, by default None
    p: float or int or None, optional
        *See* :paramref:`.CompressibleFlow.p`, by default None
    rho: float or int or None, optional
        *See* :paramref:`.CompressibleFlow.rho`, by default None
    T: float or int or None, optional
        *See* :paramref:`.CompressibleFlow.T, by default None`
    p_t: float or int or None, optional
        *See* :paramref:`.CompressibleFlow.p_t`, by default None
    rho_t: float or int or None, optional
        *See* :paramref:`.CompressibleFlow.rho_t`, by default None
    T_t: float or int or None, optional
        *See* :paramref:`.CompressibleFlow.T_t`, by default None

    Raises
    ------
    ValueError
        Insufficient inputs provided to calculate sonic velocity.

    Returns
    -------

    a : float
        Speed of sound, :math:`a`, [:math:`length/time`]
    """
    if T:
        a = (gamma*R*T)**.5

    elif p and rho:
        a = (gamma*p/rho)**0.5

    elif u == 0 and T_t:
        a = (gamma*R*T_t)**.5
            
    elif u == 0 and p_t and rho_t:
        a = (gamma*p_t/rho_t)**0.5

    else:
        raise ValueError('Lacking sufficient definition of properties required to calculate sonic velocity.')

    print(f'Speed of sound = {a}')
    return a

def get_fluid_velocity(gamma=1.4,R=287,M=None,a=None,p=None,rho=None,T=None):
    """Calculates fluid velocity, :math:`u`, based on Mach number and 
    thermodynamic properties.

    This function can solve for the fluid velocity in several different ways.

    Given :math:`M` and :math:`a`:

        * :math:`u = M*a`
    
    Given :math:`M` and :math:`T`:

        * :math:`u = M * \\sqrt{\\gamma*R*T}`

    Given :math:`M`, :math:`P`, and :math:`\\rho`:

        * :math:`u = M * \\sqrt{\\gamma*p/\\rho}`

    If :math:`M = 0`, the flow is stagnant and :math:`u = 0`.

    Parameters
    ----------
    gamma : float, optional
        *See* :paramref:`.CompressibleFlow.gamma`, by default 1.4 for 
        calorically perfect air
    R : int or float
        *See* :paramref:`.CompressibleFlow.r`, by default 287 J/kg*K
    M : float or int or None, optional
        *See* :paramref:`.CompressibleFlow.M`, by default None
    a :float or int or None, optional
        *See* :paramref:`.CompressibleFlow.a`, by default None
    p: float or int or None, optional
        *See* :paramref:`.CompressibleFlow.p`, by default None
    rho: float or int or None, optional
        *See* :paramref:`.CompressibleFlow.rho`, by default None
    T: float or int or None, optional
        *See* :paramref:`.CompressibleFlow.T`, by default None

    Returns
    -------
    u : float
        Fluid velocity, :math:`u`, [:math:`length/time`]
    """
    if M and a:
        u = M*a

    elif M and not a and T:
        u = M*(gamma*R*T)**0.5

    elif M and not T and p and rho:
        u = M*(gamma*p/rho)**0.5

    elif M == 0:
        u = 0

    print(f'Fluid velocity = {u}')
    return u

def get_mach_number(gamma=1.4,R=287,a=None,u=None,p=None,rho=None,T=None,p_t_ratio=None,rho_t_ratio=None,T_t_ratio=None):
    """Calculates fluid Mach number, :math:`M`.

    This function can solve for the fluid Mach number in several different ways.

    Given :math:`u` and :math:`a`:

        * :math:`M = u/a`
    
    Given :math:`M` and :math:`T`:

        * :math:`M = u / \\sqrt{\\gamma*R*T}`

    Given :math:`M`, :math:`p`, and :math:`\\rho`:

        * :math:`M = u / \\sqrt{\\gamma*p/\\rho}`

    Given :math:`p_t/p`:

        * :math:`M = \\sqrt{((\\frac{p_t}{p})^{\\frac{\\gamma-1}{\\gamma}} - 1) * \\frac{2}{\\gamma-1}}`

    Given :math:`\\rho_t/\\rho`:

        * :math:`M = \\sqrt{((\\frac{\\rho_t}{\\rho})^{\\gamma-1} - 1) * \\frac{2}{\\gamma-1}}`

    Given :math:`T_t/T`:

        * :math:`M = \\sqrt{((\\frac{T_t}{T}) - 1) * \\frac{2}{\\gamma-1}}`

    Parameters
    ----------
    gamma : float, optional
        *See* :paramref:`.CompressibleFlow.gamma`, by default 1.4 for 
        calorically perfect air
    R : int or float
        *See* :paramref:`.CompressibleFlow.r`, by default 287 J/kg*K 
    a :float or int or None, optional
        *See* :paramref:`.CompressibleFlow.a`, by default None
    u: float or int or None, optional
        *See* :paramref:`.CompressibleFlow.u`, by default None
    p: float or int or None, optional
        *See* :paramref:`.CompressibleFlow.p`, by default None
    rho: float or int or None, optional
        *See* :paramref:`.CompressibleFlow.rho`, by default None
    T: float or int or None, optional
        *See* :paramref:`.CompressibleFlow.T`, by default None
    p_t_ratio : float or int or None, optional
        *See* :paramref:`.CompressibleFlow.p_t_ratio`, by default None
    rho_t_ratio : float or int or None, optional
        *See* :paramref:`.CompressibleFlow.rho_t_ratio`, by default None
    T_t_ratio : float or int or None, optional
        *See* :paramref:`.CompressibleFlow.T_t_ratio`, by default None

    Returns
    -------
    M : float or int
        Mach number, :math:`u/a`
    """
    if u and a:
        M = u/a

    elif u and not a and T:
        M = u/(gamma*R*T)**0.5

    elif u and not a and not T and p and rho:
        M = u/(gamma*p/rho)**0.5

    elif p_t_ratio:
        M = ((p_t_ratio**((gamma-1)/gamma)-1)*2/(gamma-1))**0.5
        pass

    elif rho_t_ratio:
        M = ((rho_t_ratio**(gamma-1)-1)*2/(gamma-1))**0.5
        pass 

    elif T_t_ratio:
        M = ((T_t_ratio-1)*2/(gamma-1))**0.5
        pass

    elif u == 0:
        M = 0

    print(f'Mach number= {M}')
    return M

def get_static_pressure(gamma=1.4,R=287,M=None,u=None,rho=None,T=None,p_t=None,rho_t=None,T_t=None,p_t_ratio=None,rho_t_ratio=None,T_t_ratio=None):
    """Calculates fluid static pressure, :math:`p`.

    This function can solve for the fluid static pressure in several different ways.

    Given :math:`\\rho` and :math:`T`:

        * :math:`p = \\rho*R*T`
    
    Given :math:`M` and :math:`p_t`:

        * :math:`p = \\frac{p_t}{(1 + \\frac{\\gamma-1}{2} * M^2)^{\\frac{\\gamma}{\\gamma-1}}}`

    Given :math:`p_t/p` and :math:`p_t`:

       * :math:`p = \\frac{p_t}{\\frac{p_t}{p}}`

    Given :math:`p_t` and :math:`T_t/T`:

        * :math:`p = \\frac{p_t}{{(\\frac{T_t}{T}})^{\\frac{\\gamma}{\\gamma-1}}}`

    Given :math:`p_t` and :math:`T_t` and :math:`T`:

        * :math:`p = \\frac{p_t}{{(\\frac{T_t}{T}})^{\\frac{\\gamma}{\\gamma-1}}}`

    Given :math:`p_t` and :math:`\\rho_t/\\rho`:

        * :math:`p = \\frac{p_t}{{(\\frac{\\rho_t}{\\rho}})^{\\gamma}}`

    Given :math:`p_t` and :math:`\\rho_t` and :math:`\\rho`:

        * :math:`p = \\frac{p_t}{{(\\frac{\\rho_t}{\\rho}})^{\\gamma}}`

    If fluid velocity or Mach number are 0 the flow is stagnant and :math:`p = p_t`.

    Parameters
    ----------
    gamma : float, optional
        *See* :paramref:`.CompressibleFlow.gamma`, by default 1.4 for 
        calorically perfect air
    R : int or float
        *See* :paramref:`.CompressibleFlow.r`, by default 287 J/kg*K
    M : float or int or None, optional
        *See* :paramref:`.CompressibleFlow.M`, by default None
    u: float or int or None, optional
        *See* :paramref:`.CompressibleFlow.u`, by default None
    rho: float or int or None, optional
        *See* :paramref:`.CompressibleFlow.rho`, by default None
    T: float or int or None, optional
        *See* :paramref:`.CompressibleFlow.T, by default None`
    p_t: float or int or None, optional
        *See* :paramref:`.CompressibleFlow.p_t`, by default None
    rho_t: float or int or None, optional
        *See* :paramref:`.CompressibleFlow.rho_t`, by default None
    T_t: float or int or None, optional
        *See* :paramref:`.CompressibleFlow.T_t`, by default None
    p_t_ratio : float or int or None, optional
        *See* :paramref:`.CompressibleFlow.p_t_ratio`, by default None
    rho_t_ratio : float or int or None, optional
        *See* :paramref:`.CompressibleFlow.rho_t_ratio`, by default None
    T_t_ratio : float or int or None, optional
        *See* :paramref:`.CompressibleFlow.T_t_ratio`, by default None

    Returns
    -------
    p : float 
        fluid static pressure, :math:`p`, [:math:`Force/Area`]
    """
    if rho and T:
        p = rho*R*T
    
    elif p_t and M:
        p = p_t/(1 + (gamma-1) / 2 * M**2)**(gamma/(gamma-1))

    elif p_t_ratio and p_t:
        p = p_t/p_t_ratio

    elif p_t and T_t_ratio:
        p = p_t/(T_t_ratio**(gamma/(gamma-1)))

    elif p_t and T and T_t:
        p = p_t/((T_t/T)**(gamma/(gamma-1)))

    elif p_t and rho_t_ratio:
        p = p_t/rho_t_ratio**gamma

    elif p_t and rho and rho_t:
        p = p_t/(rho_t/rho)**gamma

    elif p_t and (u == 0):
        p = p_t

    print(f'Static pressure is {p}')
    return p

def get_total_pressure(gamma=1.4,R=287,M=None,u=None,p=None,rho=None,T=None,rho_t=None,T_t=None,p_t_ratio=None,rho_t_ratio=None,T_t_ratio=None):
    """Calculates fluid total pressure, :math:`p_t`.

    This function can solve for the fluid total pressure in several different ways.

    Given :math:`\\rho_t` and :math:`T_t`:

        * :math:`p_t = \\rho*R_t*T_t`
    
    Given :math:`M` and :math:`p`:

        * :math:`p_t = p*(1 + \\frac{\\gamma-1}{2} * M^2)^{\\frac{\\gamma}{\\gamma-1}}`

    Given :math:`p_t/p` and :math:`p`:

       * :math:`p_t = p*\\frac{p_t}{p}`

    Given :math:`p` and :math:`T_t/T`:

        * :math:`p_t = p*{(\\frac{T_t}{T}})^{\\frac{\\gamma}{\\gamma-1}}`

    Given :math:`p` and :math:`T_t` and :math:`T`:

        * :math:`p_t = p*{(\\frac{T_t}{T}})^{\\frac{\\gamma}{\\gamma-1}}`

    Given :math:`p` and :math:`\\rho_t/\\rho`:

        * :math:`p_t = p*{(\\frac{\\rho_t}{\\rho}})^{\\gamma}`

    Given :math:`p` and :math:`\\rho_t` and :math:`\\rho`:

        * :math:`p_t = p*{(\\frac{\\rho_t}{\\rho}})^{\\gamma}`

    If fluid velocity or Mach number are 0 the flow is stagnant and :math:`p_t = p`.

    Parameters
    ----------
    gamma : float, optional
        *See* :paramref:`.CompressibleFlow.gamma`, by default 1.4 for calorically perfect air
    R : int or float
        *See* :paramref:`.CompressibleFlow.r`, by default 287 J/kg*K
    M : float or int or None, optional
        *See* :paramref:`.CompressibleFlow.M`, by default None
    u: float or int or None, optional
        *See* :paramref:`.CompressibleFlow.u`, by default None
    p : _type_, optional
        *See* :paramref:`.CompressibleFlow.p`, by default None
    rho: float or int or None, optional
        *See* :paramref:`.CompressibleFlow.rho`, by default None
    T: float or int or None, optional
        *See* :paramref:`.CompressibleFlow.T, by default None`
    rho_t: float or int or None, optional
        *See* :paramref:`.CompressibleFlow.rho_t`, by default None
    T_t: float or int or None, optional
        *See* :paramref:`.CompressibleFlow.T_t`, by default None
    p_t_ratio : float or int or None, optional
        *See* :paramref:`.CompressibleFlow.p_t_ratio`, by default None
    rho_t_ratio : float or int or None, optional
        *See* :paramref:`.CompressibleFlow.rho_t_ratio`, by default None
    T_t_ratio : float or int or None, optional
        *See* :paramref:`.CompressibleFlow.T_t_ratio`, by default None

    Returns
    -------
    p_t : float
        fluid total pressure, :math:`p_t`, [:math:`Force/Area`]
    """
    if rho_t and T_t:
        p_t = rho_t*R*T_t

    elif p and M:
        p_t = (1 + (gamma-1) / 2 * M**2)**(gamma/(gamma-1))*p

    elif p_t_ratio and p:
        p_t = p_t_ratio*p

    elif p and T_t_ratio:
        p_t = p*(T_t_ratio**(gamma/(gamma-1)))

    elif p and T and T_t:
        p_t = p*((T_t/T)**(gamma/(gamma-1)))

    elif p and rho_t_ratio:
        p_t = p*rho_t_ratio**gamma

    elif p and rho and rho_t:
        p_t = p*(rho_t/rho)**gamma     

    elif p and (u == 0):
        p_t = p
    
    print(f'Total pressure = {p_t}')
    return p_t

def get_pressure_ratio(gamma=1.4,R=287,M=None,u=None,p=None,rho=None,T=None,p_t=None,rho_t=None,T_t=None,rho_t_ratio=None,T_t_ratio=None):
    """Calculates fluid total-to-static pressure ratio, :math:`p_t/p`.

    This function can solve for the fluid total-to-static pressure ratio in several different ways.
    
    Given :math:`M`:

        * :math:`p_t/p = (1 + \\frac{\\gamma-1}{2} * M^2)^{\\frac{\\gamma}{\\gamma-1}}`

    Given :math:`p_t` and :math:`p`:

       * :math:`p_t/p = p_t/p`

    Given :math:`\\rho_t/\\rho`:

        * :math:`p_t/p = (\\frac{\\rho_t}{\\rho})^{\\gamma}`

    Given :math:`\\rho_t` and :math:`\\rho`:

        * :math:`p_t/p = (\\frac{\\rho_t}{\\rho})^{\\gamma}`

    Given :math:`T_t/T`:

        * :math:`p_t/p = (\\frac{T_t}{T})^{\\frac{\\gamma}{\\gamma-1}}`

    Given :math:`T_t` and :math:`T`:

        * :math:`p_t/p = (\\frac{T_t}{T})^{\\frac{\\gamma}{\\gamma-1}}`

    If fluid velocity or Mach number are 0 the flow is stagnant and :math:`p_t/p = 1`.

    Parameters
    ----------
    gamma : float, optional
        *See* :paramref:`.CompressibleFlow.gamma`, by default 1.4 for 
        calorically perfect air
    R : int or float
        *See* :paramref:`.CompressibleFlow.r`, by default 287 J/kg*K
    M : float or int or None, optional
        *See* :paramref:`.CompressibleFlow.M`, by default None
    u: float or int or None, optional
        *See* :paramref:`.CompressibleFlow.u`, by default None
    p : _type_, optional
        *See* :paramref:`.CompressibleFlow.p`, by default None
    rho: float or int or None, optional
        *See* :paramref:`.CompressibleFlow.rho`, by default None
    T: float or int or None, optional
        *See* :paramref:`.CompressibleFlow.T, by default None`
    p_t: float or int or None, optional
        *See* :paramref:`.CompressibleFlow.p_t`, by default None
    rho_t: float or int or None, optional
        *See* :paramref:`.CompressibleFlow.rho_t`, by default None
    T_t: float or int or None, optional
        *See* :paramref:`.CompressibleFlow.T_t`, by default None
    rho_t_ratio : float or int or None, optional
        *See* :paramref:`.CompressibleFlow.rho_t_ratio`, by default None
    T_t_ratio : float or int or None, optional
        *See* :paramref:`.CompressibleFlow.T_t_ratio`, by default None

    Returns
    -------
    p_t_ratio : float
        total-to-static pressure ratio, :math:`p_t/p`
    """
    if M:
        p_t_ratio = (1 + (gamma-1) / 2 * M**2)**(gamma/(gamma-1))

    elif p and p_t:
        p_t_ratio = p_t/p

    elif rho_t_ratio:
        p_t_ratio = rho_t_ratio**(gamma)

    elif T_t_ratio:
        p_t_ratio = T_t_ratio**(gamma/(gamma-1))

    elif u == 0:
        p_t_ratio = 1

    print(f'p_t/p = {p_t_ratio}')
    return p_t_ratio

def get_static_temperature(gamma=1.4,R=287,M=None,u=None,p=None,rho=None,p_t=None,rho_t=None,T_t=None,p_t_ratio=None,rho_t_ratio=None,T_t_ratio=None):
    """Calculates fluid static temperature, :math:`T`.

    This function can solve for the fluid static temperature in several different ways.

    Given :math:`\\rho` and :math:`p`:

        * :math:`T = p/(\\rho*R)`
    
    Given :math:`M` and :math:`T_t`:

        * :math:`T = \\frac{T_t}{(1 + \\frac{\\gamma-1}{2} * M^2)}`

    Given :math:`T_t/T` and :math:`T_t`:

       * :math:`T = \\frac{T_t}{\\frac{T_t}{T}}`

    Given :math:`T_t` and :math:`p_t/p`:

        * :math:`T = \\frac{T_t}{{(\\frac{p_t}{p}})^{\\frac{\\gamma-1}{\\gamma}}}`

    Given :math:`T_t` and :math:`p_t` and :math:`p`:

        * :math:`T = \\frac{T_t}{{(\\frac{p_t}{p}})^{\\frac{\\gamma-1}{\\gamma}}}`

    Given :math:`T_t` and :math:`\\rho_t/\\rho`:

        * :math:`T = \\frac{T_t}{{(\\frac{\\rho_t}{\\rho}})^{\\gamma-1}}`

    Given :math:`T_t` and :math:`\\rho_t` and :math:`\\rho`:

        * :math:`T = \\frac{T_t}{{(\\frac{\\rho_t}{\\rho}})^{\\gamma-1}}`

    If fluid velocity or Mach number are 0 the flow is stagnant and :math:`T = T_t`.

    Parameters
    ----------
    gamma : float, optional
        *See* :paramref:`.CompressibleFlow.gamma`, by default 1.4 for 
        calorically perfect air
    R : int or float
        *See* :paramref:`.CompressibleFlow.r`, by default 287 J/kg*K
    M : float or int or None, optional
        *See* :paramref:`.CompressibleFlow.M`, by default None
    u: float or int or None, optional
        *See* :paramref:`.CompressibleFlow.u`, by default None
    p : _type_, optional
        *See* :paramref:`.CompressibleFlow.p`, by default None
    rho: float or int or None, optional
        *See* :paramref:`.CompressibleFlow.rho`, by default None
    p_t: float or int or None, optional
        *See* :paramref:`.CompressibleFlow.p_t`, by default None   
    rho_t: float or int or None, optional
        *See* :paramref:`.CompressibleFlow.rho_t`, by default None
    T_t: float or int or None, optional
        *See* :paramref:`.CompressibleFlow.T_t`, by default None
    p_t_ratio : float or int or None, optional
        *See* :paramref:`.CompressibleFlow.p_t_ratio`, by default None
    rho_t_ratio : float or int or None, optional
        *See* :paramref:`.CompressibleFlow.rho_t_ratio`, by default None
    T_t_ratio : float or int or None, optional
        *See* :paramref:`.CompressibleFlow.T_t_ratio`, by default None

    Returns
    -------
    T : float
        fluid static temperature, :math:`T`, [:math:`degrees`].
    """
    if p and rho:
        T = p/(rho*R)

    elif T_t and M:
        T = T_t/((1 + (gamma-1) / 2 * M**2))

    elif T_t_ratio and T_t:
        T = T_t/T_t_ratio

    elif p_t_ratio and T_t:
        T = T_t/p_t_ratio**((gamma-1)/gamma)

    elif p and p_t and T_t:
        T = T_t/(p_t/p)**((gamma-1)/gamma)

    elif rho_t_ratio and T_t:
        T = T_t/rho_t_ratio**((gamma-1))

    elif rho and rho_t and T_t:
        T = T_t/(rho_t/rho)**((gamma-1))

    elif T_t and (u == 0):
        T = T_t

    print(f'Static temperature = {T}')
    return T

def get_total_temperature(gamma=1.4,R=287,M=None,u=None,p=None,rho=None,T=None,p_t=None,rho_t=None,p_t_ratio=None,rho_t_ratio=None,T_t_ratio=None):
    """Calculates fluid total temperature, :math:`T_t`

    This function can solve for the fluid total temperature in several different ways.

    Given :math:`\\rho_t` and :math:`p_t`:

        * :math:`T_t = p_t/(\\rho_t*R)`
    
    Given :math:`M` and :math:`T`:

        * :math:`T_t = T*(1 + \\frac{\\gamma-1}{2} * M^2)`

    Given :math:`T_t/T` and :math:`T`:

       * :math:`T_t = T*\\frac{T_t}{T}`

    Given :math:`T` and :math:`p_t/p`:

        * :math:`T_t = T*(\\frac{p_t}{p})^{\\frac{\\gamma-1}{\\gamma}}`

    Given :math:`T` and :math:`p_t` and :math:`p`:

        * :math:`T_t = T*(\\frac{p_t}{p})^{\\frac{\\gamma-1}{\\gamma}}`

    Given :math:`T` and :math:`\\rho_t/\\rho`:

        * :math:`T_t = T*(\\frac{\\rho_t}{\\rho})^{\\gamma-1}`

    Given :math:`T` and :math:`\\rho_t` and :math:`\\rho`:

        * :math:`T_t = T*(\\frac{\\rho_t}{\\rho})^{\\gamma-1}`

    If fluid velocity or Mach number are 0 the flow is stagnant and :math:`T_t = T`.

    Parameters
    ----------
    gamma : float, optional
        *See* :paramref:`.CompressibleFlow.gamma`, by default 1.4 for 
        calorically perfect air
    R : int or float
        *See* :paramref:`.CompressibleFlow.r`, by default 287 J/kg*K
    M : float or int or None, optional
        *See* :paramref:`.CompressibleFlow.M`, by default None
    u: float or int or None, optional
        *See* :paramref:`.CompressibleFlow.u`, by default None
    p : _type_, optional
        *See* :paramref:`.CompressibleFlow.p`, by default None
    rho: float or int or None, optional
        *See* :paramref:`.CompressibleFlow.rho`, by default None
    T: float or int or None, optional
        *See* :paramref:`.CompressibleFlow.T, by default None`
    p_t: float or int or None, optional
        *See* :paramref:`.CompressibleFlow.p_t`, by default None   
    rho_t: float or int or None, optional
        *See* :paramref:`.CompressibleFlow.rho_t`, by default None
    p_t_ratio : float or int or None, optional
        *See* :paramref:`.CompressibleFlow.p_t_ratio`, by default None
    rho_t_ratio : float or int or None, optional
        *See* :paramref:`.CompressibleFlow.rho_t_ratio`, by default None
    T_t_ratio : float or int or None, optional
        *See* :paramref:`.CompressibleFlow.T_t_ratio`, by default None

    Returns
    -------
    T_t : float
        fluid total temperature, :math:`T_t`, [:math:`degrees`].
    """
    if p_t and rho_t:
        T_t = p_t/(rho_t*R)

    elif T and M:
        T_t = T*((1 + (gamma-1) / 2 * M**2))

    elif T and T_t:
        T_t = T*T_t_ratio

    elif p_t_ratio and T:
        T_t = T*p_t_ratio**((gamma-1)/gamma)

    elif p and p_t and T:
        T_t = T*(p_t/p)**((gamma-1)/gamma)            

    elif rho_t_ratio and T:
        T_t = T*rho_t_ratio**((gamma-1))

    elif rho and rho_t and T:
        T_t = T*(rho_t/rho_t_ratio)**((gamma-1))

    elif T and (u == 0):
        T_t = T

    print(f'Total temperature = {T_t}')
    return T_t

def get_temperature_ratio(gamma=1.4,R=287,M=None,u=None,p=None,rho=None,T=None,p_t=None,rho_t=None,T_t=None,p_t_ratio=None,rho_t_ratio=None):
    """Calculates fluid total-to-static temperature ratio, :math:`T_t/T`.

    This function can solve for the fluid total-to-static temperature ratio in several different ways.
    
    Given :math:`M`:

        * :math:`T_t/T = (1 + \\frac{\\gamma-1}{2} * M^2)`

    Given :math:`T_t` and :math:`T`:

       * :math:`T_t/T = T_t/T`

    Given :math:`\\rho_t/\\rho`:

        * :math:`T_t/T = (\\frac{\\rho_t}{\\rho})^{\\gamma-1}`

    Given :math:`\\rho_t` and :math:`\\rho`:

        * :math:`T_t/T = (\\frac{\\rho_t}{\\rho})^{\\gamma-1}`

    Given :math:`p_t/p`:

        * :math:`T_t/T = (\\frac{p_t}{p})^{\\frac{\\gamma-1}{\\gamma}}`

    Given :math:`p_t` and :math:`p`:

        * :math:`T_t/T =(\\frac{p_t}{p})^{\\frac{\\gamma-1}{\\gamma}}`

    If fluid velocity or Mach number are 0 the flow is stagnant and :math:`T_t/T = 1`.

    Parameters
    ----------
    gamma : float, optional
        *See* :paramref:`.CompressibleFlow.gamma`, by default 1.4 for 
        calorically perfect air
    R : int or float
        *See* :paramref:`.CompressibleFlow.r`, by default 287 J/kg*K
    M : float or int or None, optional
        *See* :paramref:`.CompressibleFlow.M`, by default None
    u: float or int or None, optional
        *See* :paramref:`.CompressibleFlow.u`, by default None
    p : _type_, optional
        *See* :paramref:`.CompressibleFlow.p`, by default None
    rho: float or int or None, optional
        *See* :paramref:`.CompressibleFlow.rho`, by default None
    T: float or int or None, optional
        *See* :paramref:`.CompressibleFlow.T, by default None`
    p_t: float or int or None, optional
        *See* :paramref:`.CompressibleFlow.p_t`, by default None   
    rho_t: float or int or None, optional
        *See* :paramref:`.CompressibleFlow.rho_t`, by default None
    T_t: float or int or None, optional
        *See* :paramref:`.CompressibleFlow.T_t`, by default None
    p_t_ratio : float or int or None, optional
        *See* :paramref:`.CompressibleFlow.p_t_ratio`, by default None
    rho_t_ratio : float or int or None, optional
        *See* :paramref:`.CompressibleFlow.rho_t_ratio`, by default None

    Returns
    -------
    T_t_ratio : float
        total-to-static temperature ratio, :math:`T_t/T`
    """
    if M:
        T_t_ratio = (1 + (gamma-1) / 2 * M**2)

    elif T and T_t:
        T_t_ratio = T_t/T

    elif rho_t_ratio:
        T_t_ratio = rho_t_ratio**(gamma-1)

    elif p_t_ratio:
        T_t_ratio = p_t_ratio**((gamma-1)/gamma)

    elif u == 0:
        T_t_ratio = 1

    print(f'T_t/T = {T_t_ratio}')
    return T_t_ratio

def get_static_density(gamma=1.4,R=287,M=None,u=None,p=None,T=None,p_t=None,rho_t=None,T_t=None,p_t_ratio=None,rho_t_ratio=None,T_t_ratio=None):
    """Calculates fluid static density, :math:`\\rho`

    This function can solve for the fluid static density in several different ways.

    Given :math:`T` and :math:`p`:

        * :math:`\\rho = p/(R*T)`
    
    Given :math:`M` and :math:`\\rho_t`:

        * :math:`\\rho = \\frac{\\rho_t}{(1 + \\frac{\\gamma-1}{2} * M^2)^{\\frac{1}{\\gamma-1}}}`

    Given :math:`\\rho_t/\\rho` and :math:`\\rho_t`:

       * :math:`\\rho = \\frac{\\rho_t}{\\frac{\\rho_t}{\\rho}}`

    Given :math:`\\rho_t` and :math:`p_t/p`:

        * :math:`\\rho = \\frac{\\rho_t}{{(\\frac{p_t}{p}})^{\\frac{1}{\\gamma}}}`

    Given :math:`\\rho_t` and :math:`p_t` and :math:`p`:

        * :math:`\\rho = \\frac{\\rho_t}{{(\\frac{p_t}{p}})^{\\frac{1}{\\gamma}}}`

    Given :math:`\\rho_t` and :math:`T_t/T`:

        * :math:`\\rho = \\frac{\\rho_t}{{(\\frac{T_t}{T}})^{\\frac{1}{\\gamma-1}}}`

    Given :math:`\\rho_t` and :math:`T_t` and :math:`T`:

        * :math:`\\rho = \\frac{\\rho_t}{{(\\frac{T_t}{T}})^{\\frac{1}{\\gamma-1}}}`

    If fluid velocity or Mach number are 0 the flow is stagnant and :math:`\\rho = \\rho_t`.

    Parameters
    ----------
    gamma : float, optional
        *See* :paramref:`.CompressibleFlow.gamma`, by default 1.4 for 
        calorically perfect air
    R : int or float
        *See* :paramref:`.CompressibleFlow.r`, by default 287 J/kg*K
    M : float or int or None, optional
        *See* :paramref:`.CompressibleFlow.M`, by default None
    u: float or int or None, optional
        *See* :paramref:`.CompressibleFlow.u`, by default None
    p : _type_, optional
        *See* :paramref:`.CompressibleFlow.p`, by default None
    T: float or int or None, optional
        *See* :paramref:`.CompressibleFlow.T, by default None`
    p_t: float or int or None, optional
        *See* :paramref:`.CompressibleFlow.p_t`, by default None   
    rho_t: float or int or None, optional
        *See* :paramref:`.CompressibleFlow.rho_t`, by default None
    T_t: float or int or None, optional
        *See* :paramref:`.CompressibleFlow.T_t`, by default None
    p_t_ratio : float or int or None, optional
        *See* :paramref:`.CompressibleFlow.p_t_ratio`, by default None
    rho_t_ratio : float or int or None, optional
        *See* :paramref:`.CompressibleFlow.rho_t_ratio`, by default None
    T_t_ratio : float or int or None, optional
        *See* :paramref:`.CompressibleFlow.T_t_ratio`, by default None

    Returns
    -------
    rho : float
        fluid static density, :math:`\\rho`, [:math:`mass/volume`]
    """
    if p and T:
        rho = p/(R*T)

    elif rho_t and M:
        rho = rho_t/((1 + (gamma-1) / 2 * M**2))**(1/(gamma-1))

    elif rho_t and rho_t_ratio:
        rho = rho_t/rho_t_ratio
    
    elif p_t_ratio and rho_t:
        rho = rho_t/p_t_ratio**(1/gamma)

    elif p and p_t and rho_t:
        rho = rho_t/(p_t/p)**(1/gamma)

    elif T_t_ratio and rho_t:
        rho = rho_t/T_t_ratio**(1/(gamma-1))

    elif T and T_t and rho_t:
        rho = rho_t/(T_t/T)**(1/(gamma-1))

    elif rho_t and (u == 0):
        rho = rho_t

    print(f'Static density = {rho}')
    return rho

def get_total_density(gamma=1.4,R=287,M=None,u=None,p=None,rho=None,T=None,p_t=None,T_t=None,p_t_ratio=None,rho_t_ratio=None,T_t_ratio=None):
    """Calculates fluid total density.

    This function can solve for the fluid total density in several different ways.

    Given :math:`T_t` and :math:`p_t`:

        * :math:`\\rho_t = p_t/(R*T_t)`
    
    Given :math:`M` and :math:`\\rho`:

        * :math:`\\rho_t = \\rho * (1 + \\frac{\\gamma-1}{2} * M^2)^{\\frac{1}{\\gamma-1}}`

    Given :math:`\\rho_t/\\rho` and :math:`\\rho`:

        * :math:`\\rho_t = \\rho * \\frac{\\rho_t}{\\rho}`

    Given :math:`\\rho` and :math:`p_t/p`:

        * :math:`\\rho_t = \\rho*(\\frac{p_t}{p})^{\\frac{1}{\\gamma}}`

    Given :math:`\\rho` and :math:`p_t` and :math:`p`:

        * :math:`\\rho_t = \\rho*(\\frac{p_t}{p})^{\\frac{1}{\\gamma}}`

    Given :math:`\\rho` and :math:`T_t/T`:

        * :math:`\\rho_t = \\rho*(\\frac{T_t}{T})^{\\frac{1}{\\gamma-1}}`

    Given :math:`\\rho` and :math:`T_t` and :math:`T`:

        * :math:`\\rho_t = \\rho*(\\frac{T_t}{T})^{\\frac{1}{\\gamma-1}}`

    If fluid velocity or Mach number are 0 the flow is stagnant and :math:`\\rho = \\rho_t`.

    Parameters
    ----------
    gamma : float, optional
        *See* :paramref:`.CompressibleFlow.gamma`, by default 1.4 for 
        calorically perfect air
    R : int or float
        *See* :paramref:`.CompressibleFlow.r`, by default 287 J/kg*K
    M : float or int or None, optional
        *See* :paramref:`.CompressibleFlow.M`, by default None
    u: float or int or None, optional
        *See* :paramref:`.CompressibleFlow.u`, by default None
    p : _type_, optional
        *See* :paramref:`.CompressibleFlow.p`, by default None
    rho: float or int or None, optional
        *See* :paramref:`.CompressibleFlow.rho`, by default None
    T: float or int or None, optional
        *See* :paramref:`.CompressibleFlow.T, by default None`
    p_t: float or int or None, optional
        *See* :paramref:`.CompressibleFlow.p_t`, by default None   
    T_t: float or int or None, optional
        *See* :paramref:`.CompressibleFlow.T_t`, by default None
    p_t_ratio : float or int or None, optional
        *See* :paramref:`.CompressibleFlow.p_t_ratio`, by default None
    rho_t_ratio : float or int or None, optional
        *See* :paramref:`.CompressibleFlow.rho_t_ratio`, by default None
    T_t_ratio : float or int or None, optional
        *See* :paramref:`.CompressibleFlow.T_t_ratio`, by default None

    Returns
    -------
    rho_t : float
        fluid total density, :math:`\\rho_t`, [:math:`mass/volume`]
    """
    if p_t and T_t:
        rho_t = p_t/(R*T_t)

    elif rho and M:
        rho_t = rho*((1 + (gamma-1) / 2 * M**2))**(1/(gamma-1))

    elif rho and rho_t_ratio:
        rho_t = rho*rho_t_ratio
    
    elif p_t_ratio and rho:
        rho_t = rho*p_t_ratio**(1/gamma)

    elif p and p_t and rho:
        rho_t = rho*(p_t/p)**(1/gamma)
    
    elif T_t_ratio and rho:
        rho_t = rho*T_t_ratio**(1/(gamma-1))

    elif T and T_t and rho:
        rho_t = rho*(T_t/T)**(1/(gamma-1))

    elif rho and (u == 0):
        rho_t = rho

    print(f'Total density = {rho_t}')
    return rho_t

def get_density_ratio(gamma=1.4,R=287,M=None,u=None,p=None,rho=None,T=None,p_t=None,rho_t=None,T_t=None,p_t_ratio=None,T_t_ratio=None):
    """Calculates fluid total-to-static density ratio, :math:`\\rho_t/\\rho`.

    This function can solve for the fluid total-to-static density ratio in several different ways.
    
    Given :math:`M`:

        * :math:`\\rho_t/\\rho = (1 + \\frac{\\gamma-1}{2} * M^2)^{\\frac{1}{\\gamma-1}}`

    Given :math:`\\rho_t` and :math:`\\rho`:

       * :math:`\\rho_t/\\rho = \\rho_t/\\rho`

    Given :math:`p_t/p`:

        * :math:`\\rho_t/\\rho = (\\frac{p_t}{p})^{\\frac{1}{\\gamma}}`

    Given :math:`p_t` and :math:`p`:

        * :math:`\\rho_t/\\rho = (\\frac{p_t}{p})^{\\frac{1}{\\gamma}}`

    Given :math:`T_t/T`:

        * :math:`\\rho_t/\\rho = (\\frac{T_t}{T})^{\\frac{1}{\\gamma-1}}`

    Given :math:`T_t` and :math:`T`:

        * :math:`\\rho_t/\\rho = (\\frac{T_t}{T})^{\\frac{1}{\\gamma-1}}`

    If fluid velocity or Mach number are 0 the flow is stagnant and :math:`\\rho_t/\\rho = 1`.

    Parameters
    ----------
    gamma : float, optional
        *See* :paramref:`.CompressibleFlow.gamma`, by default 1.4 for 
        calorically perfect air
    R : int or float
        *See* :paramref:`.CompressibleFlow.r`, by default 287 J/kg*K
    M : float or int or None, optional
        *See* :paramref:`.CompressibleFlow.M`, by default None
    u: float or int or None, optional
        *See* :paramref:`.CompressibleFlow.u`, by default None
    p : _type_, optional
        *See* :paramref:`.CompressibleFlow.p`, by default None
    rho: float or int or None, optional
        *See* :paramref:`.CompressibleFlow.rho`, by default None
    T: float or int or None, optional
        *See* :paramref:`.CompressibleFlow.T, by default None`
    p_t: float or int or None, optional
        *See* :paramref:`.CompressibleFlow.p_t`, by default None   
    rho_t: float or int or None, optional
        *See* :paramref:`.CompressibleFlow.rho_t`, by default None
    T_t: float or int or None, optional
        *See* :paramref:`.CompressibleFlow.T_t`, by default None
    p_t_ratio : float or int or None, optional
        *See* :paramref:`.CompressibleFlow.p_t_ratio`, by default None
    T_t_ratio : float or int or None, optional
        *See* :paramref:`.CompressibleFlow.T_t_ratio`, by default None

    Returns
    -------
    rho_t_ratio : float
        total-to-static temperature ratio, :math:`\\rho_t/\\rho`
    """
    if M:
        rho_t_ratio = (1 + (gamma-1) / 2 * M**2)**(1/(gamma-1))

    elif rho and rho_t:
        rho_t_ratio = rho_t/rho

    elif T_t_ratio:
        rho_t_ratio = T_t_ratio**(1/(gamma-1))

    elif p_t_ratio:
        rho_t_ratio = p_t_ratio**(1/gamma)

    elif u == 0:
        rho_t_ratio = 1

    print(f'rho_t/rho = {rho_t_ratio}')
    return rho_t_ratio

def get_sonic_ratios(gamma=1.4,R=287):
    """Calculates sonic to total ratios for a given :math:`\\gamma` and R.

    This function uses the isentropic ratios with M=1 to 
    output the thermodynamic ratios of sonic condition
    to total condition for :math:`p, \\rho, T`. Default values 
    are for air.

        * :math:`p^*/p_t = (\\frac{2}{\\gamma+1})^{\\frac{\\gamma}{\\gamma-1}}`

        * :math:`\\rho^*/\\rho_t= (\\frac{2}{\\gamma+1})^{\\frac{1}{\\gamma-1}}`

        * :math:`T^*/T_t = (1+\\frac{\\gamma-1}{2})^{-1}`

    Parameters
    ----------
    gamma : float, optional
        *See* :paramref:`.CompressibleFlow.gamma`, by default 1.4 for 
        calorically perfect air
    R : int or float
        *See* :paramref:`.CompressibleFlow.r`, by default 287 J/kg*K

    Returns
    -------
    p_star_p_t : float
        ratio of sonic to total pressure
    rho_star_rho_t : float
        ratio of sonic to total density
    T_star_T_t : float
        ratio of sonic to total temperature

    """

    p_star_p_t = 1/get_pressure_ratio(gamma=gamma,R=R,M=1)
    rho_star_rho_t = 1/get_density_ratio(gamma=gamma,R=R,M=1)
    T_star_T_t = 1/get_temperature_ratio(gamma=gamma,R=R,M=1)
    print(f'p*/p_t = {p_star_p_t}')
    print(f'rho*/rho_t = {rho_star_rho_t}')
    print(f'T*/T_t = {T_star_T_t}')
    return p_star_p_t, rho_star_rho_t, T_star_T_t

def get_sonic_area_ratio(gamma=1.4,M=None):
    """Calculates the sonic area ratio, :math:`A/A^*`

    This function calculates the sonic area ratio, :math:`A/A^*`,
    for a given :math:`\\gamma` and :math:`M`.

    :math:`\\frac{A}{A^*} = \\frac{1}{M} \\left[\\frac{2}{\\gamma + 1} \\left(1 + \\frac{\\gamma - 1}{2}M^2\\right) \\right] ^ { \\frac{\\gamma + 1}{2 (\\gamma - 1)}}`

    Parameters
    ----------
    gamma : float, optional
        *See* :paramref:`.CompressibleFlow.gamma`, by default 1.4 for 
        calorically perfect air
    M : float or int or None, optional
        *See* :paramref:`.CompressibleFlow.M`, by default None

    Returns
    -------
    A_A_star : float
        Sonic area ratio, :math:`A/A^*`
    """
    A_A_star = 1 / M * (2 / (gamma + 1) * (1 + (gamma - 1) / 2 * M**2) ) ** ( (gamma + 1) / (2 * (gamma - 1) ) )
    print(f'The sonic area ratio, A/A = {A_A_star}')
    return A_A_star

def get_mass_flux(gamma=1.4,R=287,M=None,p_t=None,T_t=None):
    """Calculates the mass flux, :math:`\\dot{m}/A`

    :math:`\\dot{m}/A = \\frac{p_t}{\\sqrt{R*T_t}} \\frac{\\sqrt{\\gamma}M}{\\left({1 + \\frac{\\gamma-1}{2}M^2}\\right)^{\\frac{\\gamma+1}{2(\\gamma-1)}}}`

    Parameters
    ----------
    gamma : float, optional
        *See* :paramref:`.CompressibleFlow.gamma`, by default 1.4 for 
        calorically perfect air
    R : int or float
        *See* :paramref:`.CompressibleFlow.r`, by default 287 J/kg*K
    M : float or int or None, optional
        *See* :paramref:`.CompressibleFlow.M`, by default None
    p_t: float or int or None, optional
        *See* :paramref:`.CompressibleFlow.p_t`, by default None   
    T_t: float or int or None, optional
        *See* :paramref:`.CompressibleFlow.T_t`, by default None

    Returns
    -------
    mass_flux : float
        Mass flux, :math:`\\dot{m}/A`
    """
    mass_flux = p_t/(R*T_t)**0.5 * (gamma)**0.5 * M / (1 + (gamma-1)/2 * M**2)**((gamma+1)/(2*(gamma-1)))
    print(f'Mass flux = {mass_flux}')
    return mass_flux

def get_choked_mass_flux(gamma=1.4,R=287,p_t=None,T_t=None):
    """Calculates choked mass flux (M=1), :math:`\\dot{m}/A^*`

    :math:`\\dot{m}/A^* = \\frac{p_t}{\\sqrt{T_t}} \\sqrt{\\frac{\\gamma}{R} \\left({\\frac{2}{\\gamma+1}}\\right)^{\\frac{\\gamma+1}{\\gamma-1}}}`

    Parameters
    ----------
    gamma : float, optional
        *See* :paramref:`.CompressibleFlow.gamma`, by default 1.4 for 
        calorically perfect air
    R : int or float
        *See* :paramref:`.CompressibleFlow.r`, by default 287 J/kg*K
    p_t: float or int or None, optional
        *See* :paramref:`.CompressibleFlow.p_t`, by default None   
    T_t: float or int or None, optional
        *See* :paramref:`.CompressibleFlow.T_t`, by default None

    Returns
    -------
    choked_mass_flux : float
        Mass flux through sonic (M=1) throat, :math:`\\dot{m}/A^*`
    """
    choked_mass_flux = get_mass_flux(gamma=gamma,R=R,M=1,p_t=p_t,T_t=T_t)
    print(f'Choked mass flux = {choked_mass_flux}')
    return choked_mass_flux

def get_mass_flow(gamma=1.4,R=287,M=None,p_t=None,T_t=None,A=None):
    """Calculates mass flow, :math:`\\dot{m}`, through a given area, :math:`A`

    :math:`\\dot{m} = A * \\frac{p_t}{\\sqrt{R*T_t}} \\frac{\\sqrt{\\gamma}M}{\\left({1 + \\frac{\\gamma-1}{2}M^2}\\right)^{\\frac{\\gamma+1}{2(\\gamma-1)}}}`

    Parameters
    ----------
    gamma : float, optional
        *See* :paramref:`.CompressibleFlow.gamma`, by default 1.4 for 
        calorically perfect air
    R : int or float
        *See* :paramref:`.CompressibleFlow.r`, by default 287 J/kg*K
    M : float or int or None, optional
        *See* :paramref:`.CompressibleFlow.M`, by default None
    p_t: float or int or None, optional
        *See* :paramref:`.CompressibleFlow.p_t`, by default None   
    T_t: float or int or None, optional
        *See* :paramref:`.CompressibleFlow.T_t`, by default None
    A : float or int or None, optional
        area associated with flow at a certain point of interest, by default None

    Returns
    -------
    mass_flow : float
        Mass flow rate, :math:`\\dot{m}`, through given area, :math:`A`
    """
    mass_flow = A*get_mass_flux(gamma=gamma,R=R,M=M,p_t=p_t,T_t=T_t)
    print(f'Mass flow = {mass_flow}')
    return mass_flow

def get_choked_mass_flow(gamma=1.4,R=287,p_t=None,T_t=None,A_throat=None):
    """Calculates mass flow, :math:`\\dot{m}_{max}`, through sonic (M=1) throat with given area, :math:`A_{throat}`.

    :math:`\\dot{m}_{max}` = A_{throat} * \\frac{p_t}{\\sqrt{T_t}} \\sqrt{\\frac{\\gamma}{R} \\left({\\frac{2}{\\gamma+1}}\\right)^{\\frac{\\gamma+1}{\\gamma-1}}}`

    Parameters
    ----------
    gamma : float, optional
        *See* :paramref:`.CompressibleFlow.gamma`, by default 1.4 for 
        calorically perfect air
    R : int or float
        *See* :paramref:`.CompressibleFlow.r`, by default 287 J/kg*K
    p_t: float or int or None, optional
        *See* :paramref:`.CompressibleFlow.p_t`, by default None   
    T_t: float or int or None, optional
        *See* :paramref:`.CompressibleFlow.T_t`, by default None
    A_throat : float or int or None, optional
        throat area associated with sonic (M=1) flow by default None

    Returns
    -------
    choked_mass_flow : float
        Mass flow, :math:`\\dot{m}_{max}`, through sonic (M=1) throat with given area, :math:`A_{throat}`.
    """
    choked_mass_flow = A_throat*get_choked_mass_flux(gamma=gamma,R=R,M=1,p_t=p_t,T_t=T_t)
    print(f'Choked mass flow = {choked_mass_flow}'), 
    return choked_mass_flow


if __name__=='__main__':
    #u = get_fluid_velocity(M=2,a=100)
    #u1 = get_fluid_velocity(M=2,T=300)
    #u2 = get_fluid_velocity(M=2,p=101300,rho=1.22)
    #m = get_mach_number(a=100,u=300)
    #m1 = get_mach_number(u=300,T=400)
    #m2 = get_mach_number(u=300,p=101300,rho=1.22)
    foo = CompressibleFlow(p_t_ratio=1.5)
    foo.isentropic_state()
    foo.report_outputs()
    bar = CompressibleFlow(M=1.2,p=120000,T=300)
    bar.isentropic_state()
    bar.report_outputs()
    #one = get_mach_number(T_t_ratio=2)
    #two = get_mach_number(p_t_ratio=2)
    #three = get_mach_number(rho_t_ratio=2)
    foo,bar,baz = get_sonic_ratios()
    print('I am done')