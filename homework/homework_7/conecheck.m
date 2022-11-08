function [value,isterminal,direction]=conecheck(theta,y,gamma)
% check cone solution for point where v_theta=0 (i.e., at cone surface)
% theta - current angle
% y - current solution vector
% gamma - ratio of specific heats

value=zeros(2,1);
isterminal=zeros(2,1);
direction=zeros(2,1);

%quit if Vr goes negative (non-physical solution)
value(1)=1.0;
if (y(1)<0.0)
    value(1)=0.0;
end
isterminal(1)=1;
direction(1)=0;

%quit if Vtheta goes positive (which occurs at the wall)
value(2)=1.0;
if (y(2)>0.0)
    value(2)=0.0;
end
isterminal(2)=1;
direction(2)=0;