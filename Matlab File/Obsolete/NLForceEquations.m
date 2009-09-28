function [Qx,Z] = NLForceEquations(P,Bi,Bx)
% calcualte the element forces to minimize the difference between Bx'*V
%
% written by T.Y. Yang 2009/08/23

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% use the MATLAB fmincon.m to determine the minimization function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%options = optimset('MaxIter',100,'Display','off','LargeScale','off');
Qx = 1; % initial trial
[Qx,Z] = fminunc(@(Qx)minfun(Qx,P,Bi,Bx),Qx);

%%%%%%%%%%%%%%%%%%%%%%%
% minimization function
%%%%%%%%%%%%%%%%%%%%%%%
function z = minfun(Qx,P,Bi,Bx)
Q = Bi*P+Bx*Qx;
V = Constutive(Q);
z = abs(sum(Bx'*V));

% constiutive relationship
function V = Constutive(Q)
Fy1 = 50;
E1 = 30;
b1 = 0.05;
Fy2 = 40;
E2 = 40;
b2 = 0.05;
Fy3 = 30;
E3 = 50;
b3 = 0.05;
if abs(Q(1)) > Fy1
    V1 = sign(Q(1))*(Fy1/E1+(abs(Q(1))-Fy1)/(E1*b1));
else
    V1 = Q(1)/E1;
end
if abs(Q(2)) > Fy2
    V2 = sign(Q(2))*(Fy2/E2+(abs(Q(2))-Fy2)/(E2*b2));
else
    V2 = Q(2)/E2;
end
if abs(Q(3)) > Fy3
    V3 = sign(Q(3))*(Fy3/E3+(abs(Q(3))-Fy3)/(E3*b3));
else
    V3 = Q(3)/E3;
end
V = [V1; V2; V3];