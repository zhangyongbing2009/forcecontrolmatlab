function Qx = NLForceDynamic(Pf,Bi,Bx,M,C,Uflast3,dt,Qx)
% calcualte the element forces to minimize the difference between Bx'*V
%
% written by T.Y. Yang 2009/08/23

% extract the displacement history
U1_1 = Uflast3(1,end);
U1_2 = Uflast3(1,end-1);
U1_3 = Uflast3(1,end-2);
U2_1 = Uflast3(2,end);
U2_2 = Uflast3(2,end-1);
U2_3 = Uflast3(2,end-2);
U3_1 = Uflast3(3,end);
U3_2 = Uflast3(3,end-1);
U3_3 = Uflast3(3,end-2);

tol = 1e-10;
residual = 1;

Q3 = Qx;
V3 = Constutive(Q3,3);
U2_0 = -V3;
U2dot = (1/2*(3*U2_0-4*U2_1+U2_2))/dt;
U2dotdot = (2*U2_0-5*U2_1+4*U2_2-U2_3)/dt^2;
Q2 = m2*U2dotdot+c2*U2dot+Q3+m2*Ag;
V2 = Constutive(Q2,2);
U1_0 = U2_0-V2;
U1dot = (1/2*(3*U1_0-4*U1_1+U1_2))/dt;
U1dotdot = (2*U1_0-5*U1_1+4*U1_2-U1_3)/dt^2;
Q1 = m1*U1dotdot+c1*U1dot+Q2+m1*Ag;
V1 = Constutive(Q1,1);
V = [V1,V2,V3];



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
function [V,DV] = Constutive(Q,n)
switch n
    case 1
        Fy = 15000;
        E = 2.8;
        b = 0.01;
    case 2
        Fy = 10000;
        E = 2;
        b = 0.05;
    case 3
        Fy = 10000;
        E = 5.6;
        b = 0.05;
end
if abs(Q) > Fy
    V = sign(Q)*(Fy/E+(abs(Q)-Fy)/(E*b));
    DV = 1/(b*E);
else
    V = Q/E;
    DV = 1/E;
end