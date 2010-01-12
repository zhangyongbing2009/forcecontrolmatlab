function [model analysis state] = FM_NLDynamicNRwFixIter(MODEL, ANALYSIS, STATE)
% Nonlinear dynamic analysis using force method with fixed number of iterations
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

% State Variables
U = STATE.U;
Udot = STATE.Udot;
Udotdot = STATE.Udotdot;
Ptp1 = STATE.Ptp1;
pr = STATE.pr;
prTrial = pr;
i = STATE.i;
   
% assemble the augmented matrices
Mb = [M; zeros(nos,ndf)];
Cb = [C; zeros(nos,ndf)];
Pb = [Ptp1; zeros(nos,1)];

for iter=1:maxIter

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
    Pr = B*prTrial;
    Prb = [Pr; Bx'*u];
    Sb = [B; Bx'*f];

    % get rhs and jacobian
    R = Mb*UdotdotTrial + Cb*UdotTrial + Prb - Pb;
    dRdQ = (c3*Mb + c2*Cb)*Bi'*f + Sb;

    % solve for force increments
    deltaQ = dRdQ\(-R);

    % substeps
    x = iter/maxIter;
    % scaleddeltaQ = x*(prTrial + deltaQ) - (x-1)*pr - prTrial;
    scaleddeltaQ = deltaQ/(maxIter-iter+1);

    % update variables
    prTrial = prTrial + scaleddeltaQ;

    % set trial forces in elements
    for j=1:numElem
        feval(Element{j},'setTrialStress',MatData(j),prTrial(j,1));
    end
end

% get displacements and flexibilities from elements
for j=1:numElem
    u(j,1) = feval(Element{j},'getStrain',MatData(j));
    f(j,j) = feval(Element{j},'getTangentF',MatData(j));
end

% update response quantities
UTrial = Bi'*u;
UdotTrial = c2*(UTrial-U) + a1*Udot + a2*Udotdot;
UdotdotTrial = c3*(UTrial-U) + a3*Udot + a4*Udotdot;

% commit the elements
for j=1:numElem
    feval(Element{j},'commitState',MatData(j));
end

state.U = UTrial;
state.Udot = UdotTrial;
state.Udotdot = UdotdotTrial;
state.pr = pr;
state.u = u;
state.Pr = Pr;
state.iter = maxIter;
state.errorNorm = errorNorm;
model.f = f;
model.K = MODEL.K;
analysis = ANALYSIS;
