function [model analysis state] = DM_NLDynamicNRLimitIncr(MODEL, ANALYSIS, STATE)

% Nonlinear dynamic analysis using displacement method
%
% Written by T.Y. Yang and Andreas Schellenberg and Hong Kim 09/14/2009

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
Udot = a1*STATE.Udot + a2*STATE.Udotdot;
Udotdot = a3*STATE.Udot + a4*STATE.Udotdot;
Ptp1 = STATE.Ptp1;
Pr = STATE.Pr;
pr = STATE.pr;
i = STATE.i;

% Newton-Raphson algorithm
% Initialize algorithm
iter = 0;
scaleddeltaU = zeros(ndf,1);

% get rhs and jacobian
R  = M*Udotdot + C*Udot + Pr - Ptp1;
dRdU = c3*M + c2*C + c1*K;

% solve for displacement increments
deltaU = dRdU\(-R);

if max(abs(deltaU)) > incrLimit
    scale = incrLimit/max(abs(deltaU));
    scaleddeltaU = scale*deltaU;
else
    scaleddeltaU = deltaU;
end

% update the error norm and iteration number
%errorNorm = norm(deltaU);
errorNorm = norm(R);
% errorNorm = 1.0;
while ((errorNorm >= tol) && (iter <= maxIter))
    %update response quantities
    U = U + c1*scaleddeltaU;
    Udot = Udot +c2*scaleddeltaU;
    Udotdot = Udotdot + c3*scaleddeltaU;
    
    % transform displacements from global to element DOF
    u = A*(U);
    
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
    K = A'*k*A;

    % get rhs and jacobian
    R  = M*Udotdot + C*Udot + Pr - Ptp1;
    dRdU = c3*M + c2*C + c1*K;

    % solve for displacement increments
    deltaU = dRdU\(-R);

    if max(abs(deltaU)) > incrLimit
        scale = incrLimit/max(abs(deltaU));
        scaleddeltaU = scale*deltaU;
    else
        scaleddeltaU = deltaU;
    end

    % update the error norm and iteration number
    %errorNorm = norm(deltaU);
    errorNorm = norm(R);
    iter = iter + 1;
end

state.U = U;
state.Udot = Udot;
state.Udotdot = Udotdot;
state.pr = pr;
state.u = u;
state.Pr = Pr;
state.iter = iter;
state.errorNorm = errorNorm;
model.K = K;
model.f = f;
analysis = ANALYSIS;

if (iter < maxIter)
    % commit the elements
    for j=1:numElem
%         feval(Element{j},'commitState',MatData(j));
    end
else
    error(['failed to converge in Newton-Raphson algorithm: Step = ',num2str(i),...
        ', errorNorm = ',num2str(errorNorm)]);
end
