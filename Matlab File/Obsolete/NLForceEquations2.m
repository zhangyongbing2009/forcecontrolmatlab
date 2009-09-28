function [Q,U,Udot,Udotdot,Z] = NLForceEquations2(Bx,Bi,U_1,Udot_1,Udotdot_1,m1,m2,c1,c2,dt,Ag,Qx)
% calcualte the element forces to minimize the difference between Bx'*V
%
% written by T.Y. Yang 2009/08/23

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% use the MATLAB fmincon.m to determine the minimization function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% extract the displacement history
% U1_1 = Uflast3(1,end);
% U1_2 = Uflast3(1,end-1);
% U1_3 = Uflast3(1,end-2);
% U2_1 = Uflast3(2,end);
% U2_2 = Uflast3(2,end-1);
% U2_3 = Uflast3(2,end-2);

U1_1 = U_1(1);
U1dot_1 = Udot_1(1);
U1dotdot_1 = Udotdot_1(1);
U2_1 = U_1(2);
U2dot_1 = Udot_1(2);
U2dotdot_1 = Udotdot_1(2);
gamma = 1/2;
beta = 1/4;

%Qx = 0; % initial trial
[Qx,Z] = fminunc(@(Qx)minfun(Qx,Bx,U_1,Udot_1,Udotdot_1,m1,m2,c1,c2,dt,Ag),Qx);
Q3 = Qx;
V3 = Constutive(Q3,3);
U2_0 = -V3;
[U2dot, U2dotdot] = NemarkIntegration(U2dotdot_1,U2dot_1,U2_1,U2_0,gamma,beta,dt);
%U2dot = (1/2*(3*U2_0-4*U2_1+U2_2))/dt;
%U2dotdot = (2*U2_0-5*U2_1+4*U2_2-U2_3)/dt^2;
Q2 = -m2*Ag-m2*U2dotdot-c2*U2dot-Q3;
V2 = Constutive(Q2,2);
U1_0 = U2_0-V2;
[U1dot, U1dotdot] = NemarkIntegration(U1dotdot_1,U1dot_1,U1_1,U1_0,gamma,beta,dt);
%U1dot = (1/2*(3*U1_0-4*U1_1+U1_2))/dt;
%U1dotdot = (2*U1_0-5*U1_1+4*U1_2-U1_3)/dt^2;
Q1 = -m1*Ag-m1*U1dotdot-c1*U1dot+Q2;
V1 = Constutive(Q1,1);
%V = [V1,V2,V3]';
Q = [Q1,Q2,Q3]';
%U = Bi'*V;
U = [U1_0,U2_0]';
Udot = [U1dot,U2dot]';
Udotdot = [U1dotdot,U2dotdot]';

%%%%%%%%%%%%%%%%%%%%%%%
% minimization function
%%%%%%%%%%%%%%%%%%%%%%%
function z = minfun(Qx,Bx,U_1,Udot_1,Udotdot_1,m1,m2,c1,c2,dt,Ag)
% extract the displacement history
% U1_1 = Uflast3(1,end);
% U1_2 = Uflast3(1,end-1);
% U1_3 = Uflast3(1,end-2);
% U2_1 = Uflast3(2,end);
% U2_2 = Uflast3(2,end-1);
% U2_3 = Uflast3(2,end-2);
U1_1 = U_1(1);
U1dot_1 = Udot_1(1);
U1dotdot_1 = Udotdot_1(1);
U2_1 = U_1(2);
U2dot_1 = Udot_1(2);
U2dotdot_1 = Udotdot_1(2);

gamma = 1/2;
beta = 1/4;

Q3 = Qx;
V3 = Constutive(Q3,3);
U2_0 = -V3;
[U2dot, U2dotdot] = NemarkIntegration(U2dotdot_1,U2dot_1,U2_1,U2_0,gamma,beta,dt);
% U2dot = (1/2*(3*U2_0-4*U2_1+U2_2))/dt;
% U2dotdot = (2*U2_0-5*U2_1+4*U2_2-U2_3)/dt^2;
%Q2 = m2*U2dotdot+c2*U2dot+Q3+m2*Ag;
Q2 = -m2*Ag-m2*U2dotdot-c2*U2dot-Q3;
V2 = Constutive(Q2,2);
U1_0 = U2_0-V2;
[U1dot, U1dotdot] = NemarkIntegration(U1dotdot_1,U1dot_1,U1_1,U1_0,gamma,beta,dt);
% U1dot = (1/2*(3*U1_0-4*U1_1+U1_2))/dt;
% U1dotdot = (2*U1_0-5*U1_1+4*U1_2-U1_3)/dt^2;
%Q1 = m1*U1dotdot+c1*U1dot+Q2+m1*Ag;
Q1 = -m1*Ag-m1*U1dotdot-c1*U1dot+Q2;
V1 = Constutive(Q1,1);
V = [V1,V2,V3]';
z = abs(sum(Bx'*V));

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