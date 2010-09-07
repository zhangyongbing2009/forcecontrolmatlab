function [model analysis state] = DM_NMDynamicExplicit(MODEL, ANALYSIS, STATE)

% Newmark Explicit Method
%
% Written by Hong Kim 08/24/2010

% Model Variables
A = MODEL.A;
B = MODEL.B;
C = MODEL.C;
K = MODEL.K;
M = MODEL.M;
f = MODEL.f;
Element = MODEL.Element;
MatData = MODEL.MatData;
numElem = MODEL.numElem;
ndf = MODEL.ndf;

% Analysis Variables
deltaT = ANALYSIS.deltaT;
beta = ANALYSIS.beta;
gamma = ANALYSIS.gamma;
c1 = ANALYSIS.c1;
c2 = ANALYSIS.c2;
c3 = ANALYSIS.c3;
a1 = ANALYSIS.a1;
a2 = ANALYSIS.a2;
a3 = ANALYSIS.a3;
a4 = ANALYSIS.a4;
maxIter = ANALYSIS.maxIter;
tol = ANALYSIS.tol;
incrLimit = ANALYSIS.incrLimit;

% State Variables
U = STATE.U;
u = STATE.u;
offsetu = STATE.offsetu;
Udot = STATE.Udot;
Udotdot = STATE.Udotdot;
Ptp1 = STATE.Ptp1;
Pr = STATE.Pr;
pr = STATE.pr;
i = STATE.i;

% Newton-Raphson algorithm
% Initialize algorithm
iter = 0;
errorNorm = 0;

% Calculate predictor values from Newmark approximation

UPred = U + deltaT*Udot + 0.5*deltaT^2*Udotdot;
UdotPred = Udot + deltaT*(1-gamma)*Udotdot;

u = A*UPred;
    
% Find the offset control displacment
uCtrl = u - offsetu;

% set trial response in elements
for j=1:numElem
    feval(Element{j},'setIncrTrialStrain',MatData(j),uCtrl(j,:));
end
% get resisting forces and stiffness from elements
for j=1:numElem
    pr(j,1) = feval(Element{j},'getStress',MatData(j));
    k(j,j)  = feval(Element{j},'getTangentK',MatData(j));
end

% transform forces and stiffness from element to global DOF
Pr = A'*pr;

% get rhs and jacobian
MHat = M + deltaT*gamma*C;
R = -C*UdotPred-Pr+Ptp1;

% solve for the acceleration
Udotdot = MHat\R;

state.U = UPred;
state.Udot = UdotPred + deltaT*gamma*Udotdot;
state.Udotdot = Udotdot;
state.pr = pr;
state.u = u;
state.Pr = Pr;
state.iter = iter;
state.errorNorm = errorNorm;
model.K = MODEL.K;
model.f = f;
analysis = ANALYSIS;

if (iter <= maxIter || errorNorm <= tol)
    % commit the elements
    for j=1:numElem
%          feval(Element{j},'commitState',MatData(j));
    end
else
    error(['failed to converge in Newton-Raphson algorithm: Step = ',num2str(i),...
        ', errorNorm = ',num2str(errorNorm)]);
end
