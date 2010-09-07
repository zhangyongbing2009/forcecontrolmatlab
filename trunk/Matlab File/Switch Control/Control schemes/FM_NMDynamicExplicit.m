function [model analysis state] = FM_NMDynamicExplicit(MODEL, ANALYSIS, STATE)

% Newmark Explicit Method using predictor force and measured displacement
% of the experimental element.  Then equilibrating the rest of the
% analytical portion of the structure enforcing compatibility to provide
% the resisting force of the structure in the time integration scheme.
%
% Written by Hong Kim 08/27/2010

% Model Variables
A = MODEL.A;
B = MODEL.B;
Bi = MODEL.Bi;
Bx = MODEL.Bx;
C = MODEL.C;
K = MODEL.K;
f = MODEL.f;
M = MODEL.M;
Element = MODEL.Element;
MatData = MODEL.MatData;
numElem = MODEL.numElem;
ndf = MODEL.ndf;
nos = MODEL.nos;

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
Udot = STATE.Udot;
Udotdot = STATE.Udotdot;
Ptp1 = STATE.Ptp1;
pr = STATE.pr;
offsetpr = STATE.offsetpr;
offsetu = STATE.offsetu;
i = STATE.i;
   
% Initialize algorithm
iter = 0;
errorNorm = 0;

% Calculate predictor values from Newmark approximation
UdotPred = Udot + deltaT*(1-gamma)*Udotdot;
deltaUPred = deltaT*Udot + 0.5*deltaT^2*Udotdot;
prPred = pr + inv(f)*A*deltaUPred;
pr(1,1) = prPred(1,1);
% Calculate the control forces
prCtrl = prPred - offsetpr;

% set predictor forces in experimental elements
for j=1:1
%     if isequal(Element{j},'Experimental')
        feval(Element{j},'setIncrTrialStress',MatData(j),prCtrl(j,:));
%     end
end
% get displacements and flexibilities from elements
for j=1:1
%     if isequal(Element{j},'Experimental')
        u(j,1) = feval(Element{j},'getStrain',MatData(j));
        f(j,j)  = feval(Element{j},'getTangentF',MatData(j));
%     end
end

% Formulate the nodal displacements from the measure deformation
U(1,1) = u(1,1);
U(2,1) = U(2,1) + deltaUPred(2,1);

u = A*U;

% Find the offset control displacment
uCtrl = u - offsetu;

% set predictor displacements in analytical elements
for j=2:numElem
    feval(Element{j},'setIncrTrialStrain',MatData(j),uCtrl(j,:));
end
% get resisting forces and stiffness from elements
for j=2:numElem
    pr(j,1) = feval(Element{j},'getStress',MatData(j));
    k(j,j)  = feval(Element{j},'getTangentK',MatData(j));
end

% transform forces from element to global DOF
Pr = B*pr;

MHat = M + deltaT*gamma*C;
R = -C*UdotPred-Pr+Ptp1;

% solve for the acceleration
Udotdot = MHat\R;

state.U = U;
state.Udot = UdotPred + deltaT*gamma*Udotdot;
state.Udotdot = Udotdot;
state.pr = pr;
state.u = u;
state.Pr = Pr;
state.iter = iter;
state.errorNorm = errorNorm;
model.f = f;
model.K = MODEL.K;
analysis = ANALYSIS;

if (iter <= maxIter || errorNorm <= tol)
    % commit the elements
    for j=1:numElem
%         feval(Element{j},'commitState',MatData(j));
    end
else
    error(['failed to converge in Newton-Raphson algorithm: Step = ',num2str(i),...
        ', errorNorm = ',num2str(errorNorm)]);
end


