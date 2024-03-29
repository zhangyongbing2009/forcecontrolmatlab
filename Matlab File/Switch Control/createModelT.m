function model = createModelT()

%%%%%%%%%%%%%%%%%%%%
% Model Generation
%%%%%%%%%%%%%%%%%%%%
% equilibrium matrices
B = [1 -1 0; 0 1 -1];
Bi = [1 1; 0 1; 0 0];
Bx = [1; 1; 1];

% B = [1 -1 0 -1; 0 1 -1 1];
% Bi = [1 1; 0 1; 0 0; 0 0];
% Bx = [1 0; 1 -1; 1 0; 0 1];

% B = [1 -1 0 1; 0 1 -1 0];
% Bi = [1 1; 0 1; 0 0; 0 0];
% Bx = [1 -1; 1 0; 1 0; 0 1];

% degree of static indeterminancy
nos = size(Bx,2);

A = B';

% element 1 properties
% Element{1} = 'Elastic';
Element{1} = 'BiLinearElastic';
% Element{1} = 'BiLinearHysteretic';
% Element{1} = 'Hardening';
% Element{1} = 'NLElastic';
% Element{1} = 'Experimental';
MatData(1).tag    = 1;
MatData(1).E      = 45;
MatData(1).fy     = 300;      % yield stress
MatData(1).b      = 100.0;     % hardening ratio
MatData(1).Hkin   = MatData(1).b/(1-MatData(1).b)*MatData(1).E;
MatData(1).Kiso   = 0.0;
MatData(1).amp   = 2;
MatData(1).ipAddr = '127.0.0.1';
MatData(1).ipPort = 8090;

% element 2 properties
Element{2} = 'Elastic';
MatData(2).tag = 1;
% MatData(2).E   = 200; %15
MatData(2).E   = 10;

% element 3 properties
Element{3} = 'Elastic';
MatData(3).tag = 3;
% MatData(3).E   = 1500;  %40
MatData(3).E   = 20;

% element 4 properties
Element{4} = 'Elastic';
MatData(4).tag = 4;
MatData(4).E   = 60;

% number of elements
numElem = size(B,2);

% initial stiffness matrix
k = zeros(numElem);
f = zeros(numElem);
for j=1:numElem
   feval(Element{j},'initialize',MatData(j));
   k(j,j) = feval(Element{j},'getInitialTangentK',MatData(j));
   f(j,j) = feval(Element{j},'getInitialTangentF',MatData(j));
end
K = A'*k*A;

% set the periods
T = [0.5; 0.25];
omega = (T./(2*pi)).^(-1);

% Calculate mass matrix give K and T
M = findMassMat(K,T);

% number of global DOF in the model
ndf = size(M,1);

% mass-proportional damping matrix% ElementFrc{1} = 'BiLinearHystereticForce';
zeta = 0.05;
alphaM = 2.0*zeta*omega(1);
C = alphaM*M;

% Rayleigh damping matrix
%zeta = [0.05, 0.05];
%coef = [1/omega(1) omega(1); 1/omega(2) omega(2)]\[zeta(1); zeta(2)]*2;
%alphaM = coef(1);
%betaKi = coef(2);
%C = alphaM*M + betaKi*K;

% Store model variables
model.A = A;
model.B = B;
model.Bi = Bi;
model.Bx = Bx;
model.C = C;
model.M = M;
model.K = K;
model.f = f;
model.Element = Element;
model.MatData = MatData;
model.numElem = numElem;
model.ndf = ndf;
model.nos = nos;
model.T = T;