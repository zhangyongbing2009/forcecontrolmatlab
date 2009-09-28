% Nonlinear dynamic analysis using force method
% Same as DynamicIntegratedForceMethod2.m, except fixed number of iteration is used 
%
% Written by T.Y. Yang and Andreas Schellenberg 09/08/2009

% clean start
clear all; close all; clc

% location to save the trial force
fileName = 'Qtrial.txt'; %[];%'Qtrial.txt';
FID = fopen(fileName,'a');

% load the ground motion
path = 'C:\Documents and Settings\Tony\Desktop\Force control\Models\ThreeSpanTruss\';
% get ground-acceleration
ag = load([path,'\elcentro.txt']);
ag = 386.1*ag;
dt = 0.02;
% load OpenSees outputs
data = load([path,'\Node_Dsp.out']);
tN = [0; data(:,1)];
dN = [0 0; data(:,[2,3])];
data = load([path,'\Node_Vel.out']);
vN = [0 0; data(:,[2,3])];
data = load([path,'\Node_Acc.out']);
aN = [0 0; data(:,[2,3])];
% get element data
data = load([path,'\Elmt_Frc.out']);
fE = [0 0 0; data(:,2:4)];
% get element data
dataGF = load([path,'\Elmt_GFrc.out']);
fGE = [0 0 0 0; dataGF(:,[4,6,8,10])];
data = load([path,'\Elmt_Defo.out']);
dE = [0 0 0; data(:,2:4)];

Q1O = fE(:,1); % tension positive
Q2O = fE(:,2); % tension positive
Q3O = fE(:,3); % tension positive
Q1GO = fGE(:,1); % tension positive
Q2GLO = fGE(:,2); % tension positive
Q2GRO = fGE(:,3); % tension positive
Q3GO = fGE(:,4); % tension positive

V1O = dE(:,1); % tension positive
V2O = dE(:,2); % tension positive
V3O = dE(:,3); % tension positive

U1O = dN(:,1);
U2O = dN(:,2);
U1dotO = vN(:,1);
U2dotO = vN(:,2);
U1dotdotO = aN(:,1);
U2dotdotO = aN(:,2);


% COnstants for the three span truss problem
M = [0.04 0; 0 0.02];
C = 1.01*M;
B = [1 -1 0; 0 1 -1];
Bi = [1 1; 0 1; 0 0];
Bx = [1; 1; 1];
%Ft = diag([1/2.8,1/2, 1/5.6]);
[V1,DV1] = Constutive(0,1,FID);
[V2,DV2] = Constutive(0,2);
[V3,DV3] = Constutive(0,3);
Ft = diag([DV1,DV2,DV3]);

Ib = [1; 1];

% Numark constants
gamma = 1/2;
beta = 1/4;
c2 = gamma/beta/dt;
c3 = 1/beta/dt/dt;
a1 = 1-gamma/beta;
a2 = dt*(1-1/2*gamma/beta);
a3 = -1/beta/dt;
a4 = 1-1/2/beta;

Tol = 1e-5;
maxIter = 10;

% intial trial of Q
Q = zeros(size(Bi,1),length(ag));
V = zeros(size(Bi,1),length(ag));
U = zeros(size(M,2),length(ag)); 
Udot = zeros(size(M,2),length(ag)); 
Udotdot = zeros(size(M,2),length(ag)); 

% Assemble the matrix
Sb = [B; Bx'*Ft];
Mb = [M; zeros(size(Bx,2),size(M,2))];
Cb = [C; zeros(size(Bx,2),size(C,2))];

% calculate the structural dynamics
for nn = 2:length(ag)
    iter = 1;
    normdQ = 1;
    Pr = B*Q(:,nn);
    Prb = [Pr; zeros(size(Bx,2),1)];
    Pb = -Mb*Ib*ag(nn);
    
    % calculate the initial trial responses
    Q(:,nn) = Q(:,nn-1);    
    U(:,nn) = U(:,nn-1);
    Udot(:,nn) = a1*Udot(:,nn-1)+a2*Udotdot(:,nn-1);
    Udotdot(:,nn) = a3*Udot(:,nn-1)+a4*Udotdot(:,nn-1);
    
    for iter=1:maxIter+1
        % update the tangent flexibility
        [V1,DV1] = Constutive(Q(1,nn),1,FID);
        [V2,DV2] = Constutive(Q(2,nn),2);
        [V3,DV3] = Constutive(Q(3,nn),3);
        V(:,nn) = [V1; V2; V3];
        Ft = diag([DV1,DV2,DV3]);
        %V(:,nn) = Ft*Q(:,nn);
        
        U(:,nn) = Bi'*V(:,nn);
        Udot(:,nn) = c2*(U(:,nn)-U(:,nn-1))+a1*Udot(:,nn-1)+a2*Udotdot(:,nn-1);
        Udotdot(:,nn) = c3*(U(:,nn)-U(:,nn-1))+a3*Udot(:,nn-1)+a4*Udotdot(:,nn-1);
        
        if iter < maxIter+1
            Pr = B*Q(:,nn);
            Prb = [Pr; zeros(size(Bx,2),1)];
            Sb = [B; Bx'*Ft];
            
            % assemble the Jacobian matrix
            R = Mb*Udotdot(:,nn)+Cb*Udot(:,nn)+Prb-Pb;
            dRdQ = c3*Mb*Bi'*Ft+c2*Cb*Bi'*Ft+Sb;
            dQ = -dRdQ\R;
            x = iter/maxIter;
            scaleddQ = x*(Q(:,nn) + dQ) - (x-1)*Q(:,nn-1) - Q(:,nn);
            
            % update variables
            Q(:,nn) = Q(:,nn) + scaleddQ;
            
            % update the tolerance and iteration number
            normdQ = norm(dQ);
        end
    end
end

% figure;
subplot(3,1,1)
plot(V(1,:),Q(1,:));
grid
subplot(3,1,2)
plot(V(2,:),Q(2,:));
grid
subplot(3,1,3)
plot(V(3,:),Q(3,:));
grid
 
% plot the figure
figure;
subplot(3,1,1)
plot(Q(1,:))
hold on
plot(Q1O,':r')
ylabel('Q1')
grid
subplot(3,1,2)
plot(Q(2,:))
hold on
plot(Q2O,':r')
ylabel('Q2')
grid
subplot(3,1,3)
plot(Q(3,:))
hold on
plot(Q3O,':r')
ylabel('Q3')
xlabel('Step')
grid

figure;
subplot(3,2,1)
plot(U(1,:))
hold on
plot(U1O,':r')
ylabel('U1')
grid
subplot(3,2,2)
plot(U(2,:))
hold on
ylabel('U2')
plot(U2O,':r')
grid
subplot(3,2,3)
plot(Udot(1,:))
hold on
plot(U1dotO,':r')
ylabel('U1dot')
grid
subplot(3,2,4)
plot(Udot(2,:))
hold on
ylabel('U2dot')
plot(U2dotO,':r')
grid
subplot(3,2,5)
plot(Udotdot(1,:))
hold on
plot(U1dotdotO,':r')
ylabel('U1dotdot')
xlabel('Step')
grid
subplot(3,2,6)
plot(Udotdot(2,:))
hold on
ylabel('U2dotdot')
xlabel('Step')
plot(U2dotdotO,':r')
grid



fclose(FID);