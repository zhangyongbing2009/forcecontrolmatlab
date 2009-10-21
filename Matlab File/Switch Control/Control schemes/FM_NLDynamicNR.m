function [model analysis state] = FM_NLDynamicNR(MODEL, ANALYSIS, STATE)
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
M = MODEL.M;
Element = MODEL.Element;
MatData = MODEL.Mat;
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

% State Variables
U = STATE.U;
Udot = STATE.Udot;
Udotdot = STATE.Udotdot;
Ptp1 = STATE.Ptp1;
pr = STATE.pr;
i = STATE.i;

% initialize element variables
% u = zeros(numElem,npts);
% pr = zeros(numElem,npts);
% uall = zeros(numElem,npts);
% prall = zeros(numElem,npts);

% assemble the augmented matrices
Mb = [M; zeros(nos,ndf)];
Cb = [C; zeros(nos,ndf)];
Pb = [Ptp1; zeros(nos,1)];

% Newton-Raphson algorithm
iter = 0;
errorNorm = 1.0;
while (errorNorm >= tol && iter <= maxIter)

    % get displacements and flexibilities from elements
    for j=1:numElem
        u(j,1) = feval(Element{j},'getStrain',MatData(j));
        f(j,j) = feval(Element{j},'getTangentF',MatData(j));
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

    % update response quantity
    pr = pr + deltaQ;

    % set trial forces in elements
    for j=1:numElem
        feval(Element{j},'setTrialStress',MatData(j),pr(j,1));
    end

    % update the error norm and iteration number
    errorNorm = norm(deltaQ);
    %errorNorm = norm(R);
    iter = iter+1;
end

state.U = UTrial;
state.Udot = UdotTrial;
state.Udotdot = UdotdotTrial;
state.pr = pr;
state.u = u;
state.Pr = Pr;
model.f = f;
model.K = MODEL.K;
analysis = ANALYSIS;

if (iter < maxIter)
    % commit the elements
    for j=1:numElem
        feval(Element{j},'commitState',MatData(j));
    end
else
    error(['failed to converge in Newton-Raphson algorithm: Step = ',num2str(i),...
        ', errorNorm = ',num2str(errorNorm)]);
end


