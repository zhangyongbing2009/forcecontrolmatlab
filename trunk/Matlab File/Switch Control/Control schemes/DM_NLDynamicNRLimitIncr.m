function [model analysis state] = DM_NLDynamicNRLimit(MODEL, ANALYSIS, STATE)

% Nonlinear dynamic analysis using displacement method
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
incrLimit = ANALYSIS.incrLimit;

% State Variables
U = STATE.U;
UIncr = zeros(MODEL.ndf,1);
Udot = STATE.Udot;
Udotdot = STATE.Udotdot;
Udot = a1*STATE.Udot + a2*STATE.Udotdot;
Udotdot = a3*STATE.Udot + a4*STATE.Udotdot;
Ptp1 = STATE.Ptp1;
i = STATE.i;

% Newton-Raphson algorithm
iter = 0;
errorNorm = 1.0;
while ((errorNorm >= tol) && (iter <= maxIter))

    % get resisting forces and stiffness from elements
    for j=1:numElem
        pr(j,1) = feval(Element{j},'getStress',MatData(j));
        k(j,j)    = feval(Element{j},'getTangentK',MatData(j));
    end

    % transform forces and stiffness from element to global DOF
    Pr = A'*pr;
    K = A'*k*A;

    % get rhs and jacobian
    R  = M*Udotdot + C*Udot + Pr - Ptp1;
    dRdU = c3*M + c2*C + c1*K;
    normK = norm(K);
    normdd = norm(c3*M + c2*C);
    
    % solve for displacement increments
    deltaU = dRdU\(-R);

    if max(abs(deltaU)) > incrLimit
        scale = incrLimit/max(abs(deltaU));
        scaleddeltaU = scale*deltaU;
    else
        scaleddeltaU = deltaU;
    end

    % update response quantities
    UIncr = UIncr + c1*scaleddeltaU;
    Udot = Udot +c2*scaleddeltaU;
    Udotdot = Udotdot + c3*scaleddeltaU;

    % transform displacements from global to element DOF
    uIncr = A*UIncr;

    % set trial response in elements
    for j=1:numElem
        feval(Element{j},'setIncrTrialStrain',MatData(j),uIncr(j,:));
    end

    % update the error norm and iteration number
    errorNorm = norm(deltaU);
    %errorNorm = norm(R);
    iter = iter + 1;
end

state.U = U + UIncr;
state.Udot = Udot;
state.Udotdot = Udotdot;
state.pr = pr;
state.u = A*(U + UIncr);
state.Pr = Pr;
state.iter = iter;
state.errorNorm = errorNorm;
state.normK = normK;
state.normdd = normdd;
state.k = k;
model.K = K;
analysis = ANALYSIS;

if (iter <= maxIter || errorNorm <= tol)
    % commit the elements
    for j=1:numElem
        feval(Element{j},'commitState',MatData(j));
    end
else
    error(['failed to converge in Newton-Raphson algorithm: Step = ',num2str(nn),...
        ', errorNorm = ',num2str(errorNorm)]);
end
