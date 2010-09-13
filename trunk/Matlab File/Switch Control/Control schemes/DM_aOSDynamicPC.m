function [model analysis state] = DM_aOSDynamicPC(MODEL, ANALYSIS, STATE)

% alpha OS
%
% Written by Hong Kim 09/08/2010

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
alpha = ANALYSIS.alpha;
beta = ANALYSIS.beta;
gamma = ANALYSIS.gamma;

maxIter = ANALYSIS.maxIter;
tol = ANALYSIS.tol;

% State Variables
U = STATE.U;
UPredPrev = STATE.UPred;
u = STATE.u;
offsetu = STATE.offsetu;
Udot = STATE.Udot;
Udotdot = STATE.Udotdot;
Ptp0 = STATE.Ptp0;
Ptp1 = STATE.Ptp1;
Pr = STATE.Pr;
PrPrev = Pr;
pr = STATE.pr;
i = STATE.i;

% Newton-Raphson algorithm
% Initialize algorithm
iter = 0;
errorNorm = 0;

% Calculate predictor values from Newmark approximation
UPred = U + deltaT*Udot + (0.5 - beta)*deltaT^2*Udotdot;

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
KHat = (1/beta/deltaT^2)*M + (alpha*gamma/beta/deltaT)*gamma*C +alpha*K;
R = -alpha*(C*(Udot+deltaT*(1-gamma)*Udotdot) + Pr - Ptp1) ...
    -(1-alpha)*(C*Udot+(PrPrev + K*U-K*UPredPrev) - Ptp0);

% solve for the acceleration
deltaU = KHat\R;

state.UPred = UPred;
state.U = UPred+deltaU;
state.Udot = gamma/beta/deltaT*deltaU + Udot + deltaT*(1-gamma)*Udotdot;
state.Udotdot = 1/beta/deltaT^2 * deltaU;
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
