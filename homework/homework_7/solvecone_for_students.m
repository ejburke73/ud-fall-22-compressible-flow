function [theta_c,M_c,sol]=solvecone(theta_s,M_inf,gamma)
% numerical solover of Taylor-Maccoll eqn. (Eqn. 10.13 in Anderson)
% theta_s - shock angle [degrees]
% M_inf - upstream Mach number
% gamma - ratio of specific heats

% convert to radians and make sure the shock angle is greater than the Mach
% angle
theta_sr=theta_s.*pi./180.0;
if (theta_sr<=asin(1./M_inf))
    theta_c=0.0; M_c=M_inf; Mn1=nan; Mn2=nan;
    return;
end

% calculate initial flow deflection (delta) just after shock...
% call on theta-beta-M function to calculated delta
delta=theta_beta_M(theta_sr,M_inf,gamma)
% calculate the normal components of the Mach number pre- and post-shock
Mn1=M_inf.*sin(theta_sr);
Mn2=sqrt((Mn1.^2+(2./(gamma-1)))./(2.*gamma./(gamma-1).*Mn1.^2-1));
% calculate M2 just behind the shock
M2=Mn2./(sin(theta_sr-delta))

% calculate the non-dimensional velocity just after the oblique shock using
% eqn. 10.16
V_prime=(1+2./((gamma-1).*M2.^2)).^(-0.5);
% calculate velocity components from V_prime in spherical coordinates
Vr_prime=V_prime.*cos(theta_sr-delta);
Vtheta_prime=-V_prime.*sin(theta_sr-delta);
% calculate initial values for derivatives from eqn. 10.14
dVr_prime=Vtheta_prime;
% assemble initial value for ODE solver
y0=[Vr_prime;dVr_prime]

% set up ODE solver
% stop procedure when v_theta=0
% see conecheck fuction
% set relative error tolerance small enough to handle low M
options=odeset('Events',@conecheck,'RelTol',1e-12);
% solve by marching solution away from shock until either 0 degrees or
% flow tangency reached as indicated by y(2)==0.
% see cone function
[sol]=ode15s(@cone,[theta_sr 1e-10],y0,options,gamma);
% check if we have a solution, as ode15s may not converge for some values.
[n,m]=size(sol.ye);
theta_c=0.0;
M_c=M_inf;
% if ODE solver worked, calculate the angle and Mach number at the cone
if (n>0 & m>0 & abs(sol.ye(2))<1e-10)
    theta_c=sol.xe.*180.0./pi;
    Vc2=sol.ye(1).^2+sol.ye(2).^2;
    M_c=((1.0./Vc2-1).*(gamma-1)./2).^-0.5;
end