% Stiffness and flexibility test
% The purpose of this study is to see the effect of flexibility and
% stiffneess matrix on the number of iterations to solve an equation
clear all; clc; close all;

M = diag([1 0.5]);
% equilibrium matrices
B = [1 -1 0; 0 1 -1];
Bi = [1 1; 0 1; 0 0];
Bx = [1; 1; 1];

% Compatibility matrix
A = B';

% Stiffness matrix
k=diag([40 15 1000]);
f=inv(k);
K=A'*k*A;

K = A'*k*A;

% % set the periods
% T = [0.5; 0.25];
% omega = (T./(2*pi)).^(-1);
% 
% % Calculate mass matrix give K and T
% M = findMassMat(K,T);

% number of global DOF in the model
ndf = size(M,1);
nos = size(Bx,2);
numElem = 3;
deltaT = 0.02;

% calculate natural frequencies and periods
lambda = eig(K,M);
omega  = sort(sqrt(lambda));
T = 2.0*pi./omega;

% Viscous damping matrix
zeta = 0.05;
alphaM = 2.0*zeta*omega(1);
C = alphaM*M;

% analysis parameters
beta  = 1/4;
gamma = 1/2;
c1 = 1.0;
c2 = gamma/(beta*deltaT);
c3 = 1.0/(beta*deltaT*deltaT);
a1 = (1.0 - gamma/beta);
a2 = (deltaT)*(1.0 - 0.5*gamma/beta);
a3 = -1.0/(beta*deltaT);
a4 = 1.0 - 0.5/beta;

% Load
Ptp1 = 386.1*M*[1; 1];

% DISPLACEMENT METHOD
% max iterations and tol
maxIter_D = 1000;
tol_D = 1.0E-8;
incrLimit_D = 5E-3;

% State Variables
U_D = zeros(ndf,1);
Udot_D = zeros(ndf,1);
Udotdot_D = zeros(ndf,1);
Udot_D = a1*Udot_D + a2*Udotdot_D;
Udotdot_D = a3*Udot_D + a4*Udotdot_D;

% Newton-Raphson algorithm
iter_D = 0;
errorNorm_D = 1.0;
scaleddeltaU = zeros(ndf,1);
while ((errorNorm_D >= tol_D) && (iter_D <= maxIter_D))
   %update response quantities
    U_D = U_D + c1*scaleddeltaU;
    Udot_D = Udot_D +c2*scaleddeltaU;
    Udotdot_D = Udotdot_D + c3*scaleddeltaU;
    
    % transform displacements from global to element DOF
    u_D = A*(U_D);
   
    % get resisting forces and stiffness from elements
    for j=1:numElem
        pr_D(j,1) = k(j,j)*u_D(j);
    end

    % transform forces and stiffness from element to global DOF
    Pr_D = A'*pr_D;
    K = A'*k*A;

    % get rhs and jacobian
    R_D(:,iter_D+1)  = M*Udotdot_D + C*Udot_D + Pr_D - Ptp1;
    cM = c3*M;  cC = c2*C;  cK = c1*K;
    dRdU = cM + cC + cK;

    % solve for displacement increments
    deltaU(:,iter_D+1) = dRdU\(-R_D(:,iter_D+1));
    deltau(:,iter_D+1) = A*deltaU(:,iter_D+1);

    if abs(deltaU(1,iter_D+1)) > incrLimit_D
%     if max(abs(deltaU(:,iter_D+1))) > incrLimit_D
%         scale = incrLimit_D/max(abs(deltaU(:,iter_D+1)));
        scale = incrLimit_D/abs(deltaU(1,iter_D+1));
        scaleddeltaU = scale*deltaU(:,iter_D+1);
    else
        scaleddeltaU = deltaU(:,iter_D+1);
    end

    % update the error norm and iteration number
%     errorNorm_D = norm(deltaU(:,iter_D+1));
    errorNorm_D = norm(R_D(:,iter_D+1));
    iter_D = iter_D + 1;
end

% FORCE METHOD
% assemble the augmented matrices
Mb = [M; zeros(nos,ndf)];
Cb = [C; zeros(nos,ndf)];
Pb = [Ptp1; zeros(nos,1)];

% zero all the state variables again
U_F = zeros(ndf,1);
Udot_F = zeros(ndf,1);
Udotdot_F = zeros(ndf,1);
pr_F = zeros(numElem,1);

% max iterations and tol
maxIter_F = 1000;
tol_F = tol_D;
incrLimit_F = incrLimit_D*k(1,1);

% Newton-Raphson algorithm
iter_F = 0;
errorNorm_F = 1.0;
scaleddeltaQ = zeros(numElem,1);
while (errorNorm_F >= tol_F && iter_F <= maxIter_F)
    % update response quantity
    pr_F = pr_F + scaleddeltaQ;

    % get displacements and flexibilities from elements
    for j=1:numElem
        u_F(j,1) = f(j,j)*pr_F(j);
    end

    % update response quantities
    UTrial = Bi'*u_F;
    UdotTrial = c2*(UTrial-U_F) + a1*Udot_F + a2*Udotdot_F;
    UdotdotTrial = c3*(UTrial-U_F) + a3*Udot_F + a4*Udotdot_F;
    
    % transform forces from element to global DOF
    Pr_F = B*pr_F;
    Prb = [Pr_F; Bx'*u_F];
    Sb = [B; Bx'*f];

    % get rhs and jacobian
    R_F(:,iter_F+1) = Mb*UdotdotTrial + Cb*UdotTrial + Prb - Pb;
    cMb=c3*Mb*Bi'*f;    cCb = c2*Cb*Bi'*f;
    dRdQ = cMb + cCb + Sb;

    % solve for force increments
    deltaQ(:,iter_F+1) = dRdQ\(-R_F(:,iter_F+1));

    % Scale increment
%     if max(abs(deltaQ(:,iter_F+1))) > incrLimit_F
    if abs(deltaQ(1,iter_F+1)) > incrLimit_F
%         scale = incrLimit_F/max(abs(deltaQ(:,iter_F+1)));
        scale = incrLimit_F/abs(deltaQ(1,iter_F+1));
        scaleddeltaQ = scale*deltaQ(:,iter_F+1);
    else
        scaleddeltaQ = deltaQ(:,iter_F+1);
    end

    % update the error norm and iteration number
%     errorNorm_F = norm(deltaQ(:,iter_F+1));
    errorNorm_F = norm(B*R_F(:,iter_F+1));
    iter_F = iter_F+1;
end
U_F = UTrial;
Udot_F = UdotTrial;
Udotdot_F = UdotdotTrial;

err = U_D-U_F

iter_D
iter_F

figure
plot((R_D(1,:)'),'-+')
hold on
plot(R_D(2,:)','r-o')
plot((R_F(1,:)'),'g-+')
plot(R_F(2,:)','m-o')
grid on
title('Residuals')

figure
plot((deltau(1,:)'),'-+')
hold on
plot(deltau(2,:)','r-o')
grid on
plot(deltau(2,:)','g-+')
title('\Delta u')
legend('u_1','u_2','u_3')

figure
plot((deltaQ(1,:)'),'-+')
hold on
plot(deltaQ(2,:)','r-o')
plot(deltaQ(3,:)','g-+')
grid on
title('\Delta Q')
legend('Q_1','Q_2','Q_3')
