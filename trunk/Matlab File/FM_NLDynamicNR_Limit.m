% Nonlinear dynamic analysis using force method
%
% Written by T.Y. Yang and Andreas Schellenberg 09/08/2009

% clean start
clear all; close all; clc;

% add the subroutine path to the folder
% addpath([pwd '/Material models']);
addpath([pwd '\Material models']);

% constants for the three span truss problem
%M = [0.04 0; 0 0.02];
M = [0.01 0; 0 0.01];
%C = 1.01*M;
B = [1 -1 0; 0 1 -1];
Bi = [1 1; 0 1; 0 0];
Bx = [1; 1; 1];

% B = [1 -1 0 -1; 0 1 -1 1];
% Bi = [1 1; 0 1; 0 0; 0 0];
% Bx = [1 0; 1 -1; 1 0; 0 1];

% B = [1 -1 0 1; 0 1 -1 0];
% Bi = [1 1; 0 1; 0 0; 0 0];
% Bx = [1 -1; 1 0; 1 0; 0 1];


% number of global DOF in the model
ndf = size(M,1);
% degree of static indeterminancy
nos = size(Bx,2);

% element 1 properties - left column
%Element{1} = 'ElasticForce';
%Element{1} = 'BiLinearElasticForce';
 %Element{1} = 'BiLinearHystericForce';
%Element{1} = 'HardeningForce';
%Element{1} = 'NLElasticForce';
Element{1} = 'ExperimentalForce';
MatData(1).tag    = 1;
MatData(1).E      = 2.65;
MatData(1).Fy     = 1.5;    % yield stress
MatData(1).b      = 0.02;    % hardening ratio
MatData(1).ipAddr = '127.0.0.1';
MatData(1).ipPort = 8090;
MatData(1).Hkin   = 0.01;


% element 2 properties
Element{2} = 'ElasticForce';
MatData(2).tag = 2;
% MatData(2).E   = 2;
MatData(2).E   = 5;

% element 3 properties
Element{3} = 'ElasticForce';
MatData(3).tag = 3;
% MatData(3).E   = 5.6;
MatData(3).E   = 20;

% element 4 properties - spring
Element{4} = 'ElasticForce';
MatData(4).tag = 4;
MatData(4).E   = 10;

% initial flexibility matrix
numElem = size(B,2);
Ft = zeros(numElem);
for j=1:numElem
   feval(Element{j},'initialize',MatData(j));
   Ft(j,j) = feval(Element{j},'getInitialTangent',MatData(j));
end

K = B*inv(Ft)*B';

% calculate natural frequencies and periods
lambda = eig(K,M);
omega  = sort(sqrt(lambda));
T = 2.0*pi./omega;
%C = 1.01*M;

% Calculate Raleigh Damping
%zeta = [0.03 0.03];
zeta = [0.025 0.025];
a_o = zeta(1) * 2 * omega(1) * omega(2) / (omega(1)+omega(2));
a_1 = zeta(1) * 2 / (omega(1)+omega(2));
C = a_o*M + a_1*K;

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load GroundMotion Data
%%%%%%%%%%%%%%%%%%%%%%%%%%
% load the ground motion
GMDir = 'D:\Force Control\Ground motions\';
% GMDir = '/Users/hongkim/Research/Force Control/forcecontrolmatlab/Ground motions/';
dt = 0.02;
SF = 1.25;
g = 386.1;
ag0 = load([GMDir 'elcentro.txt']);
t0 = 0:length(ag0)-1;
t0 = dt*t0;
tEnd = t0(end);
ag0 = SF*g*ag0;

% change to analysis deltaT
deltaT = 0.02;
t = deltaT*(0:floor(tEnd/deltaT))';
ag = interp1(t0,ag0,t);
b = [1; 1];
npts = 303; %length(ag);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Newmark Transient Analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

% max iterations and Tol
maxIter = 200;
Tol = 1e-2;
incrLimit = 0.025;
counter = 1;

% initialize global response variables
Q = zeros(numElem,npts);
Qm = zeros(1,npts);
errorNorms = zeros(1,npts);
V = zeros(numElem,npts);
U = zeros(ndf,npts);
Udot = zeros(ndf,npts);
Udotdot = zeros(ndf,npts);
numIter = zeros(1,npts);

% assemble the augmented matrices
Sb = [B; Bx'*Ft];
Mb = [M; zeros(nos,ndf)];
Cb = [C; zeros(nos,ndf)];

% calculations for each time step, nn
for nn = 2:npts
   iter = 1;
   errorNorm = 1;
   Pr = B*Q(:,nn);
   Prb = [Pr; zeros(nos,1)];
   Pb = -Mb*b*ag(nn);
   
   % calculate the initial trial responses
   Q(:,nn) = Q(:,nn-1);
   U(:,nn) = U(:,nn-1);
   Udot(:,nn) = a1*Udot(:,nn-1) + a2*Udotdot(:,nn-1);
   Udotdot(:,nn) = a3*Udot(:,nn-1) + a4*Udotdot(:,nn-1);
   
   while (errorNorm >= Tol && iter <= maxIter)
      % set trial forces in elements
      for j=1:numElem
         feval(Element{j},'setTrialStress',MatData(j),Q(j,nn));
      end
      
      % get displacements and flexibilities from elements
      for j=1:numElem
         V(j,nn) = feval(Element{j},'getStrain',MatData(j));
         Ft(j,j) = feval(Element{j},'getTangent',MatData(j));
      end
      
      Qall(:,counter) = Q(:,nn);
      Vall(:,counter) = V(:,nn);
      
      U(:,nn) = Bi'*V(:,nn);
      Udot(:,nn) = c2*(U(:,nn)-U(:,nn-1)) + a1*Udot(:,nn-1) + a2*Udotdot(:,nn-1);
      Udotdot(:,nn) = c3*(U(:,nn)-U(:,nn-1)) + a3*Udot(:,nn-1) + a4*Udotdot(:,nn-1);
      
      Pr = B*Q(:,nn);
%       Prb = [Pr; zeros(nos,1)];
      Prb = [Pr; Bx'*V(:,nn)];
      Sb = [B; Bx'*Ft];
      
      % assemble the Jacobian matrix
      R = Mb*Udotdot(:,nn) + Cb*Udot(:,nn) + Prb - Pb;
      dRdQ = c3*Mb*Bi'*Ft + c2*Cb*Bi'*Ft + Sb;
      dQ = -dRdQ\R;
      
      if max(abs(dQ)) > incrLimit
         scale = incrLimit/max(abs(dQ));
         dQ = scale*dQ;
      end
      
      
      % update variables
      Q(:,nn) = Q(:,nn) + dQ;
      
      % update the tolerance and iteration number
      % get energy increment
      %errorNorm = norm(dQ);
      errorNorm = norm(R);
      iter = iter+1;
      counter = counter+1;
      
      clc
      nn
      iter
   end
   
   numIter(nn) = iter;
   errorNorms(nn) = errorNorm;
   
   if (iter < maxIter)
      % commit the elements
      for j=1:numElem
         feval(Element{j},'commitState',MatData(j));
      end
   else
      error(['failed to converge in Newton-Raphson algorithm: Step = ',num2str(nn),...
         ', errorNorm = ',num2str(errorNorm)]);
   end
   
   % set trial forces in elements
   for j=1:numElem
       feval(Element{j},'setTrialStress',MatData(j),Q(j,nn));
   end
   
   % get displacements and flexibilities from elements
   for j=1:numElem
      if isequal(Element{j},'ExperimentalForce')
         [V(j,nn),Qm(j,nn)] = feval(Element{j},'getStrain',MatData(j)); 
      else
         V(j,nn) = feval(Element{j},'getStrain',MatData(j)); 
      end            
       Ft(j,j) = feval(Element{j},'getTangent',MatData(j));
   end
   
   U(:,nn) = Bi'*V(:,nn);
   Udot(:,nn) = c2*(U(:,nn)-U(:,nn-1)) + a1*Udot(:,nn-1) + a2*Udotdot(:,nn-1);
   Udotdot(:,nn) = c3*(U(:,nn)-U(:,nn-1)) + a3*Udot(:,nn-1) + a4*Udotdot(:,nn-1);
   
   % commit the elements
   for j=1:numElem
       feval(Element{j},'commitState',MatData(j));
   end
   
end

% disconnect from experimental sites
for j=1:numElem
   if isequal(Element{j},'ExperimentalForce')
      feval(Element{j},'disconnect',MatData(j));
   end
end

%%%%%%%%%%%%%%%%%%
% plot the figures
%%%%%%%%%%%%%%%%%%
t = t(1:nn);

% plot the element hysteresis
figure;
for j=1:numElem
    subplot(numElem,1,j);
    plot(V(j,:),Q(j,:));
    xlabel(['V' num2str(j)]);
    ylabel(['Q' num2str(j)]);
    grid
end

figure;
for j=1:numElem
    subplot(numElem,1,j);
    plot(Vall(j,:),Qall(j,:));
    xlabel(['Vall' num2str(j)]);
    ylabel(['Qall' num2str(j)]);
    grid
end

% plot the element force history
figure;
for j=1:numElem
    subplot(numElem,1,j);
    plot(t,Q(j,:));
    xlabel('Time [sec]')
    ylabel(['Q' num2str(j)]);
    grid
end

% plot the element displacement history
figure;
for j=1:numElem
    subplot(numElem,1,j);
    plot(t,V(j,:));
    xlabel('Time [sec]')
    ylabel(['V' num2str(j)]);
    grid
end

% plot the Node response history
figure;
for j=1:ndf
    subplot(3,ndf,j);
    plot(t,U(j,:));
    xlabel('Time [sec]')
    ylabel(['U' num2str(j)]);
    grid    
    subplot(3,ndf,j+ndf);
    plot(t,Udot(j,:));
    xlabel('Time [sec]')
    ylabel(['Udot' num2str(j)]);
    grid
    subplot(3,ndf,j+2*ndf);
    plot(t,Udot(j,:));
    xlabel('Time [sec]')
    ylabel(['Udot' num2str(j)]);
    grid
end

% error
figure;
plot(t,errorNorms)
ylabel('Error')
xlabel('Time [sec]')
grid

figure;
plot(t,numIter)
ylabel('# of Iter/dt')
xlabel('Time [sec]')
grid

% ground motion
figure;
plot(t,ag0(1:length(t)))
ylabel('Ag')
xlabel('Time [sec]')
grid

% remove the subroutine path to the folder
% rmpath([pwd '/Material models']);
rmpath([pwd '\Material models']);