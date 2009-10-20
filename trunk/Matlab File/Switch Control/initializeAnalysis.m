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

% max iterations and tol
maxIter = 100;
tol = 1.0E-6;

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

% max iterations and tol
analysis.maxIter = maxIter;
analysis.tol = tol;