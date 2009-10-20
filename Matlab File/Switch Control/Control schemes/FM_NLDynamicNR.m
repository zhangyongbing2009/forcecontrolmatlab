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

numElem = size(B,2);


% initialize element variables
% u = zeros(numElem,npts);
% pr = zeros(numElem,npts);
% uall = zeros(numElem,npts);
% prall = zeros(numElem,npts);

% assemble the augmented matrices
Mb = [M; zeros(nos,ndf)];
Cb = [C; zeros(nos,ndf)];
   
% get new response quantity
pr(:,i+1) = pr(:,i);

% get applied loads
Pb = [Ptp1; zeros(nos,1)];

% Newton-Raphson algorithm
iter = 0;
errorNorm = 1.0;
while (errorNorm >= tol && iter <= maxIter)

    % get displacements and flexibilities from elements
    for j=1:numElem
        u(j,1) = feval(Element{j},'getStrain',MatData(j));
        f(j,j)  = feval(Element{j},'getTangentF',MatData(j));
    end

    % update response quantities
    U(:,i+1) = Bi'*u(:,i+1);
    Udot(:,i+1) = c2*(U(:,i+1)-U(:,i)) + a1*Udot(:,i) + a2*Udotdot(:,i);
    Udotdot(:,i+1) = c3*(U(:,i+1)-U(:,i)) + a3*Udot(:,i) + a4*Udotdot(:,i);

    % transform forces from element to global DOF
    Pr(:,i+1) = B*pr(:,i+1);
    Prb = [Pr(:,i+1); Bx'*u(:,i+1)];
    Sb = [B; Bx'*f];

    % get rhs and jacobian
    R(:,count) = Mb*Udotdot(:,i+1) + Cb*Udot(:,i+1) + Prb - Pb;
    dRdQ = (c3*Mb + c2*Cb)*Bi'*f + Sb;

    % solve for force increments
    deltaQ = dRdQ\(-R(:,count));

    % update response quantity
    pr(:,i+1) = pr(:,i+1) + deltaQ;
    prall(:,count+1) = pr(:,i+1);
    % set trial forces in elements
    for j=1:numElem
        feval(Element{j},'setTrialStress',MatData(j),pr(j,i+1));
    end

    % update the error norm and iteration number
    errorNorm = norm(deltaQ);
    %errorNorm = norm(R);
    iter = iter+1;
end

if (iter < maxIter)
    % commit the elements
    for j=1:numElem
        feval(Element{j},'commitState',MatData(j));
    end
else
    error(['failed to converge in Newton-Raphson algorithm: Step = ',num2str(i),...
        ', errorNorm = ',num2str(errorNorm)]);
end


