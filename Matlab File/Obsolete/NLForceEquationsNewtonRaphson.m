function [Q,U,Udot,Udotdot] = NLForceEquationsNewtonRaphson(U_1,Udot_1,Udotdot_1,m1,m2,c1,c2,dt,Ag,Qx,Tol,iterMax)
% calcualte the element forces to minimize the difference between Bx'*V
%
% written by T.Y. Yang 2009/08/23

% constants
gamma = 1/2;
beta = 1/4;
b2 = gamma/beta/dt;
b3 = 1/beta/dt/dt;

% past values
U1_1 = U_1(1);
U1dot_1 = Udot_1(1);
U1dotdot_1 = Udotdot_1(1);
U2_1 = U_1(2);
U2dot_1 = Udot_1(2);
U2dotdot_1 = Udotdot_1(2);

% initial trial
norm_dQx = 1;
iter = 1;

while norm_dQx > Tol && iter <= iterMax
    Q3 = Qx;
    [V3,F3] = Constutive(Q3,3);
    U2_0 = -V3;
    [U2dot, U2dotdot] = NemarkIntegration(U2dotdot_1,U2dot_1,U2_1,U2_0,gamma,beta,dt);
    
    Q2 = -m2*Ag-m2*U2dotdot-c2*U2dot-Q3;
    [V2,F2] = Constutive(Q2,2);
    U1_0 = U2_0-V2;
    [U1dot, U1dotdot] = NemarkIntegration(U1dotdot_1,U1dot_1,U1_1,U1_0,gamma,beta,dt);
    
    Q1 = -m1*Ag-m1*U1dotdot-c1*U1dot+Q2;
    [V1,F1] = Constutive(Q1,1);
    
    G = V1+V2+V3;
    DG = F1*(-1 + m2*b3*F3 + c2*b2*F3 + m1*b3*F3 - m1*b3*F2 + m1*b3*F2*m2*b3*F3...
        + m1*b3*F2*c2*b2*F3 + c1*b2*F3 - c1*b2*F2 + c1*b2*F2*m2*b3*F3 + c1*b2*F2*c2*b2*F3) ...
        + F2*(-1 + m2*b3*F3 + c2*b2*F3) + F3;
    dQx = -G/DG;
    Qx = Qx+dQx;
    
    norm_dQx = norm(dQx);
    
    iter = iter + 1;
end

if iter > iterMax
    error(['iter reached iterMax = ' num2str(iterMax) ', norm_dQx = ' num2str(norm_dQx)]);
end

% update the Q, U, Udot, Udotdot
Q3 = Qx;
V3 = Constutive(Q3,3);
U2_0 = -V3;
[U2dot, U2dotdot] = NemarkIntegration(U2dotdot_1,U2dot_1,U2_1,U2_0,gamma,beta,dt);
Q2 = -m2*Ag-m2*U2dotdot-c2*U2dot-Q3;
V2 = Constutive(Q2,2);
U1_0 = U2_0-V2;
[U1dot, U1dotdot] = NemarkIntegration(U1dotdot_1,U1dot_1,U1_1,U1_0,gamma,beta,dt);
Q1 = -m1*Ag-m1*U1dotdot-c1*U1dot+Q2;
Q = [Q1,Q2,Q3]';
U = [U1_0,U2_0]';
Udot = [U1dot,U2dot]';
Udotdot = [U1dotdot,U2dotdot]';