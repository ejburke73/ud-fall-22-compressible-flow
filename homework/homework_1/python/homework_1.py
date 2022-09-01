#!/usr/bin/env python3

rho = 1.23
p = 1.01e5
du_u = 0.01
gamma = 1.4

# Problem 1c
u_c = 63
drho_rho_c = -rho/(gamma*p)*u_c**2*du_u
print(f'The fractional density change is {round(drho_rho_c*100,3)} %')

# Problem 1d
u_d = 980
drho_rho_d = -rho/(gamma*p)*u_d**2*du_u
print(f'The fractional density change is {round(drho_rho_d*100,3)} %')

magnitude = drho_rho_d/drho_rho_c
print(f'The fractional density change in part (d) is {round(magnitude,0)} times that of part (c)')