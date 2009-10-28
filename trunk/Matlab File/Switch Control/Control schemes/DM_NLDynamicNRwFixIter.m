function [model analysis state] = DM_NLDynamicNRwFixIter(MODEL, ANALYSIS, STATE)
% Nonlinear dynamic analysis using displacement method with fixed number of iterations
%
% Written by T.Y. Yang and Andreas Schellenberg and Hong Kim 09/14/2009

% Model Variables
A = MODEL.A;
B = MODEL.B;
C = MODEL.C;
K = MODEL.K;
M = MODEL.M;
Element = MODEL.Element;
MatData = MODEL.MatData;
numElem = MODEL.numElem;

% Analysis Variables
c1 = ANALYSIS.c1;
c2 = ANALYSIS.c2;
c3 = ANALYSIS.c3;
a1 = ANALYSIS.a1;
a2 = ANALYSIS.a2;
a3 = ANALYSIS.a3;
a4 = ANALYSIS.a4;
maxIter = ANALYSIS.maxIter;
tol = ANALYSIS.tol;

% State Variables
U = STATE.U;
U = UTrial;
Udot = STATE.Udot;
Udotdot = STATE.Udotdot;
Udot = a1*STATE.Udot + a2*STATE.Udotdot;
Udotdot = a3*STATE.Udot + a4*STATE.Udotdot;
Ptp1 = STATE.Ptp1;
i = STATE.i;
   
for iter=1:maxIter

    % get resisting forces and stiffness from elements
    for j=1:numElem
        pr(j,1) = feval(Element{j},'getStress',MatData(j));
        k(j,j)    = feval(Element{j},'getTangent',MatData(j));
    end

    % transform forces and stiffness from element to global DOF
    Pr = A'*pr;
    K = A'*k*A;

    % get rhs and jacobian
    R  = M*Udotdot + C*Udot + Pr - Ptp1;
    dRdU = c3*M + c2*C + c1*K;

    % solve for displacement increments
    deltaU = dRdU\(-R);

    % sub-steps
    x = iter/maxIter;
    scaleddeltaU = x*(UTrial + deltaU) - (x-1)*U - UTrial;
    %scaleddeltaU = deltaU/(maxIter-iter+1);

    % update response quantities
    UTrial = UTrial + c1*scaleddeltaU;
    Udot = Udot +c2*scaleddeltaU;
    Udotdot = Udotdot + c3*scaleddeltaU;

    % transform displacements from global to element DOF
    u = A*UTrial;

    % set trial response in elements
    for j=1:numElem
        feval(Element{j},'setTrialStrain',MatData(j),u(j,1));
    end
end

% get resisting forces and stiffness from elements
for j=1:numElem
    pr(j,1) = feval(Element{j},'getStress',MatData(j));
    k(j,j)    = feval(Element{j},'getTangent',MatData(j));
end

% transform forces and stiffness from element to global DOF
Pr = A'*pr;
K = A'*k*A;

% get rhs and jacobian
R  = M*Udotdot + C*Udot + Pr - Ptp1;
dRdU = c3*M + c2*C + c1*K;

% solve for displacement increments
deltaU = dRdU\(-R);
errorNorms(i) = norm(deltaU);

% update response quantities
UTrial = UTrial + c1*scaleddeltaU;
Udot = Udot +c2*scaleddeltaU;
Udotdot = Udotdot + c3*scaleddeltaU;

% commit the elements
for j=1:numElem
    feval(Element{j},'commitState',MatData(j));
end

state.U = UTrial;
state.Udot = Udot;
state.Udotdot = Udotdot;
state.pr = pr;
state.u = u;
state.Pr = Pr;
model.K = K;
analysis = ANALYSIS;