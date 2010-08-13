function [model analysis state] = FM_NLDynamicNRLimitIncr(MODEL, ANALYSIS, STATE)
% Nonlinear dynamic analysis using force method
%
% Written by T.Y. Yang and Andreas Schellenberg 09/08/2009

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
i = STATE.i;

% assemble the augmented matrices
Mb = [M; zeros(nos,ndf)];
Cb = [C; zeros(nos,ndf)];
Pb = [Ptp1; zeros(nos,1)];
   
% Newton-Raphson algorithm
% Initialize algorithm
iter = 0;
scaleddeltaQ = zeros(numElem,1);
% update response quantities
UTrial = Bi'*u;
UdotTrial = c2*(UTrial-U) + a1*Udot + a2*Udotdot;
UdotdotTrial = c3*(UTrial-U) + a3*Udot + a4*Udotdot;

% transform forces from element to global DOF
Pr = B*pr;
Prb = [Pr; Bx'*u];
Sb = [B; Bx'*f];

% get rhs and jacobian
R = Mb*UdotdotTrial + Cb*UdotTrial + Prb - Pb;
dRdQ = (c3*Mb + c2*Cb)*Bi'*f + Sb;

% solve for force increments
deltaQ = dRdQ\(-R);

% Scale increment
if max(abs(deltaQ)) > incrLimit
    scale = incrLimit/max(abs(deltaQ));
    scaleddeltaQ = scale*deltaQ;
else
    scaleddeltaQ = deltaQ;
end

% update the error norm and iteration number
% errorNorm = norm(deltaQ);
errorNorm = norm(B*R);
while (errorNorm >= tol && iter <= maxIter)
    % update response quantity
    pr = pr + scaleddeltaQ;
    % find the offset control force
    prCtrl = pr - offsetpr;
    % set trial forces in elements
    for j=1:numElem
        feval(Element{j},'setIncrTrialStress',MatData(j),prCtrl(j,:));
    end
    % get displacements and flexibilities from elements
    for j=1:numElem
        u(j,1) = feval(Element{j},'getStrain',MatData(j));
        f(j,j)  = feval(Element{j},'getTangentF',MatData(j));
    end

    % update response quantities
    UTrial = Bi'*u;
    UdotTrial = c2*(UTrial-U) + a1*Udot + a2*Udotdot;
    UdotdotTrial = c3*(UTrial-U) + a3*Udot + a4*Udotdot;
    
    % transform forces from element to global DOF
    Pr = B*pr;
    Prb = [Pr; Bx'*u];
    Sb = [B; Bx'*f];

    % get rhs and jacobian
    R = Mb*UdotdotTrial + Cb*UdotTrial + Prb - Pb;
    dRdQ = (c3*Mb + c2*Cb)*Bi'*f + Sb;

    % solve for force increments
    deltaQ = dRdQ\(-R);

    % Scale increment
    if max(abs(deltaQ)) > incrLimit
        scale = incrLimit/max(abs(deltaQ));
        scaleddeltaQ = scale*deltaQ;
    else
        scaleddeltaQ = deltaQ;
    end

    % update the error norm and iteration number
    %errorNorm = norm(deltaQ);
    errorNorm = norm(B*R);
    iter = iter+1;
end

state.U = UTrial;
state.Udot = UdotTrial;
state.Udotdot = UdotdotTrial;
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


