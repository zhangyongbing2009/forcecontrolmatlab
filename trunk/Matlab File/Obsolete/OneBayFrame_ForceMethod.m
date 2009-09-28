
close all;
clear all;
clc;

%%%%%%%%%%%%%%%%%%%%
% Model Generation
%%%%%%%%%%%%%%%%%%%%
% number of global DOF in the model
ndf = 2;
% degree of static indeterminancy
nos = 1;

% define mass matrix
M = zeros(ndf,ndf);
M(1,1) = 0.04;   % mass for DOF 1 (kips/g)
M(2,2) = 0.02;   % mass for DOF 2 (kips/g)
Mb = [M; zeros(nos,ndf)];

% element 1 properties - left column
Element{1} = 'ElasticInv';
MatData(1).tag    = 1;
MatData(1).F      = 1/2.8;
MatData(1).ipAddr = '127.0.0.1';
MatData(1).ipPort = 8090;

% element 2 properties - right column
Element{2} = 'ElasticInv';
MatData(2).tag = 2;
MatData(2).F   = 1/5.6;

% element 3 properties - spring
Element{3} = 'ElasticInv';
MatData(3).tag = 3;
MatData(3).F   = 1/2.0;

% transformation matrices
B = [1  0 -1;
     0 -1  1];
Bi = [1 1;
      0 0;
      0 1];
Bx = [1;
      1;
      1];

% initial stiffness matrix
numElem = length(Element);
f = zeros(numElem);
% for j=1:numElem
%    feval(Element{j},'initialize',MatData(j));
%    k(j,j) = feval(Element{j},'getInitialTangent',MatData(j));
% end
% K = A'*k*A;
% 
% % calculate natural frequencies and periods
% lambda = eig(K,M);
% omega  = sort(sqrt(lambda));
% T = 2.0*pi./omega;

% mass-proportional damping matrix
zeta = 0.05;
%alphaM = 2.0*zeta*omega(1);
alphaM = 1.01;
C = alphaM*M;
Cb = [C; zeros(nos,ndf)];

% Rayleigh damping matrix
%zeta = [0.05, 0.05];
%coef = [1/omega(1) omega(1); 1/omega(2) omega(2)]\[zeta(1); zeta(2)]*2;
%alphaM = coef(1);
%betaKi = coef(2);
%C = alphaM*M + betaKi*K;


%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load GroundMotion Data
%%%%%%%%%%%%%%%%%%%%%%%%%%
dt = 0.02;
g = 386.1;
[t0 ag0] = ReadWriteTHFile('readTHF','elcentro.txt');
t0 = dt*t0;
tEnd = t0(end);
ag0 = g*ag0;

% change to analysis deltaT
deltaT = 0.02;
t = deltaT*(0:floor(tEnd/deltaT))';
ag = interp1(t0,ag0,t);
b = [1; 1];
npts = length(ag);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Newmark Transient Analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% analysis parameters
beta  = 0.25;
gamma = 0.50;
tol = 1.0E-5;
maxIter = 25;

% initialize global response variables
U = zeros(ndf,npts);
Udot = zeros(ndf,npts);
Udotdot = zeros(ndf,npts);
Pr = zeros(ndf,npts);

% initialize element variables
u = zeros(numElem,npts);
pr = zeros(numElem,npts);

% calculations for each time step, i
for i=1:npts-1
   
   % set the constants
   c1 = 1.0;
   c2 = gamma/(beta*deltaT);
   c3 = 1.0/(beta*deltaT*deltaT);
   a1 = (1.0 - gamma/beta);
   a2 = (deltaT)*(1.0 - 0.5*gamma/beta);
   a3 = -1.0/(beta*deltaT);
   a4 = 1.0 - 0.5/beta;
   
   % get new response quantities
   pr(:,i+1) = pr(:,i);
   U(:,i+1) = U(:,i);
   Udot(:,i+1) = a1*Udot(:,i) + a2*Udotdot(:,i);
   Udotdot(:,i+1) = a3*Udot(:,i) + a4*Udotdot(:,i);
   
   % get applied loads
   Ptb = -Mb*b*ag(i+1);
   
   % Newton-Raphson algorithm
   iter = 0;
   normQIncr = 1.0;
   while ((normQIncr >= tol) && (iter <= maxIter))
            
      % set trial forces in elements
      for j=1:numElem
         feval(Element{j},'setTrialStress',MatData(j),pr(j,i+1));
      end
      
      % get displacements and flexibilities from elements
      for j=1:numElem
         u(j,i+1) = feval(Element{j},'getStrain',MatData(j));
         f(j,j)   = feval(Element{j},'getTangent',MatData(j));
      end
      
      % transform displacements and flexibilities from element to global DOF
      Pr(:,i+1) = B*pr(:,i+1);
      Prb = [Pr(:,i+1); zeros(nos,1)];
      Sb = [B; Bx'*f];
      
      % get rhs and jacobian
      rhs = Mb*Udotdot(:,i+1) + Cb*Udot(:,i+1) + Prb - Ptb;
      jac = c3*Mb*Bi'*f + c2*Cb*Bi'*f + Sb;
      
      % solve for redundant force increments
      deltaQ = jac\(-rhs);
      pr(:,i+1) = pr(:,i+1) + deltaQ;
      
      % set trial forces in elements
      for j=1:numElem
         feval(Element{j},'setTrialStress',MatData(j),pr(j,i+1));
      end
      
      % get displacements and flexibilities from elements
      for j=1:numElem
         u(j,i+1) = feval(Element{j},'getStrain',MatData(j));
         f(j,j)   = feval(Element{j},'getTangent',MatData(j));
      end
      
      % update response quantities
      U(:,i+1) = Bi'*u(:,i+1);
      Udot(:,i+1) = c2*(U(:,i+1) - U(:,i)) + a1*Udot(:,i) + a2*Udotdot(:,i);
      Udotdot(:,i+1) = c3*(U(:,i+1) - U(:,i)) + a3*Udot(:,i) + a4*Udotdot(:,i);
      
      % get norm of Qx increment
      normQIncr = norm(deltaQ);
      iter = iter + 1;
   end
   
   if (iter < maxIter)
      % commit the elements
      for j=1:numElem
         feval(Element{j},'commitState',MatData(j));
      end
   else
      error('failed to converge in Newton-Raphson algorithm');
   end
end

% disconnect from experimental sites
for j=1:numElem
   if isequal(Element{j},'Experimental')
      feval(Element{j},'disconnect',MatData(j));
   end
end


%%%%%%%%%%%%%%%%%%%
% Post-Processing
%%%%%%%%%%%%%%%%%%%
% plot displacement vs. time
figure(1);
plot(t,U(1,:),'b-');
hold('on');
plot(t,U(2,:),'r-');
grid('on');
ylabel('Displacement [in.]');
xlabel('Time [sec]');
legend('DOF 1','DOF 2',1);

% plot element force vs. element deformation
figure(2);
plot(u(1,:),pr(1,:),'b-');
hold('on');
plot(u(2,:),pr(2,:),'r-');
grid('on');
ylabel('Resisting Force [kips]');
xlabel('Displacement [in.]');
legend('Element 1','Element 2',2);
