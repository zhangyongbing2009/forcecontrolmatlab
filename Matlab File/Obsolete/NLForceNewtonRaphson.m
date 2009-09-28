function [Q,Uf_0] = NLForceNewtonRaphson(Pf,Bi,Bx,M,C,Uflast3,dt,Qx)
% calcualte the element forces to minimize the difference between Bx'*V
%
% written by T.Y. Yang 2009/08/23

% extract the displacement history
U_1 = Uflast3(:,end);
U_2 = Uflast3(:,end-1);
U_3 = Uflast3(:,end-2);

tol = 1e-10;
residual = 1;

P = Pf;
Q = Bi*P+Bx*Qx;
[V,DV] = Constutive(Q);

while residual > tol
    FF = Bx'*DV*Bx;
    dQx = -FF\(Bx'*V);
    Qx = Qx+dQx;
    Q = Bi*P+Bx*Qx;
    [V,DV] = Constutive(Q);
    residual = abs(Bx'*V);
end

Uf_0 = Bi'*V;
Ufdot = (1/2*(3*Uf_0-4*U_1+U_2))/dt;
Ufdotdot = (2*Uf_0-5*U_1+4*U_2-U_3)/dt^2;
P = Pf - M*Ufdotdot - C*Ufdot;

while residual > tol
    FF = Bx'*DV*Bx;
    dQx = -FF\(Bx'*V);
    Qx = Qx+dQx;
    Q = Bi*P+Bx*Qx;
    [V,DV] = Constutive(Q);
    residual = abs(Bx'*V);
    Uf_0 = Bi'*V;
    Ufdot = (1/2*(3*Uf_0-4*U_1+U_2))/dt;
    Ufdotdot = (2*Uf_0-5*U_1+4*U_2-U_3)/dt^2;
    P = Pf - M*Ufdotdot - C*Ufdot;
end

Q = Bi*P+Bx*Qx;
V = Constutive(Q);
Uf_0 = Bi'*V;

% constiutive relationship
function [V,DV] = Constutive(Q)
Fy1 = 15000;
E1 = 2.8;
b1 = 0.01;
Fy2 = 10000;
E2 = 2;
b2 = 0.05;
Fy3 = 10000;
E3 = 5.6;
b3 = 0.05;
if abs(Q(1)) > Fy1
    V1 = sign(Q(1))*(Fy1/E1+(abs(Q(1))-Fy1)/(E1*b1));
    DV1 = 1/(b1*E1);
else
    V1 = Q(1)/E1;
    DV1 = 1/E1;
end
if abs(Q(2)) > Fy2
    V2 = sign(Q(2))*(Fy2/E2+(abs(Q(2))-Fy2)/(E2*b2));
    DV2 = 1/(b2*E2);
else
    V2 = Q(2)/E2;
    DV2 = 1/E2;
end
if abs(Q(3)) > Fy3
    V3 = sign(Q(3))*(Fy3/E3+(abs(Q(3))-Fy3)/(E3*b3));
    DV3 = 1/(b3*E3);
else
    V3 = Q(3)/E3;
    DV3 = 1/E3;
end
V = [V1; V2; V3];
DV = diag([DV1,DV2,DV3]);