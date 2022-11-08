function [dy]=cone(theta,y,gamma)
% y is a vector containing vr, vr’
% governing equations are continuity, irrotationality, & Euler’s equation.
dy=zeros(2,1);

dy(1)=y(2);
dy(2)=(y(2).^2.*y(1)-(gamma-1)./2.*(1-y(1).^2-y(2).^2).*(2.*y(1)+y(2).*cot(theta)))./((gamma-1)./2.*(1-y(1).^2-y(2).^2)-y(2).^2);
