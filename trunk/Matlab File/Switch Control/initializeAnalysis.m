function analysis = initializeAnalysis(deltaT)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Newmark Transient Analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% analysis parameters
beta  = 1/4;
gamma = 1/2;
c1 = 1.0;
c2 = gamma/(beta*deltaT);
c3 = 1.0/(beta*deltaT*deltaT);
a1 = (1.0 - gamma/beta);
a2 = (deltaT)*(1.0 - 0.5*gamma/beta);
a3 = -1.0/(beta*deltaT);
a4 = 1.0 - 0.5/beta;

% Analysis Method
% schemeDisp = 'DM_NLDynamicNR';
% schemeDisp = 'DM_NLDynamicNRLimit';
% schemeDisp = 'DM_NLDynamicNRwFixIter';
schemeDisp = 'DM_NLDynamicNRLimitIncr';
% schemeForce = 'FM_NLDynamicNR';
% schemeForce = 'FM_NLDynamicNRLimit';
% schemeForce = 'FM_NLDynamicNRwFixIter';
schemeForce = 'FM_NLDynamicNRLimitIncr';

% Switch Method
% schemeSwitch = 'dispCtrlOnly';
% schemeSwitch = 'forceCtrlOnly';
schemeSwitch = 'simpleYield';
% schemeSwitch = 'secantUpdate';

% Switch Parameters
Kd = 1.1;
Kf = 1.5;
Rd = 1.2;
Rf = 1.4;


plotFlag = 'r-';
% plotFlag = 'b-';

% max iterations and tol
maxIter = 500;
tol = 1.0E-3;
incrLimit = 5E-2;

% Store analysis variables
analysis.beta  = beta;
analysis.gamma = gamma;
analysis.c1 = c1;
analysis.c2 = c2;
analysis.c3 = c3;
analysis.a1 = a1;
analysis.a2 = a2;
analysis.a3 = a3;
analysis.a4 = a4;

% Store switching parameters
analysis.Kd  = Kd;
analysis.Kf = Kf;
analysis.Rd = Rd;
analysis.Rf = Rf;

% max iterations and tol
analysis.schemeDisp = schemeDisp;
analysis.schemeForce = schemeForce;
analysis.schemeSwitch = schemeSwitch;
analysis.plotFlag = plotFlag;
analysis.maxIter = maxIter;
analysis.tol = tol;
analysis.incrLimit = incrLimit;