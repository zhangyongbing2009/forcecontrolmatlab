function [C,Z] = Coeff_Minimization0(stress1,strain1,E,stress2,strain2,minC,maxC)
% calcualte the coefficient by fitting the RambergOsgood coeffiecnets
% 
% written by T.Y. Yang on 2009/09/09

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% use the MATLAB fmincon.m to determine the minimization function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%[C,Z] = fminunc(@(C)minfun(C,stress1,strain1,E,stress2,strain2),C);
options = optimset('MaxIter',100,'Display','off','LargeScale','off');
C = [1,1]; % initial trial
[C,Z] = fmincon(@(C)minfun(C,stress1,strain1,E,stress2,strain2), C, [],[],[],[],[],[],@(C)confun(C,minC,maxC),options);

%%%%%%%%%%%%%%%%%%%%%%%
% minimization function
%%%%%%%%%%%%%%%%%%%%%%%
function z = minfun(C,stress1,strain1,E,stress2,strain2)
A = C(1);
R = C(2);
z = abs(strain1 - stress1/E*(1 + A*(stress1/E)^(R-1))) + ...
    abs(strain2 - stress2/E*(1 + A*(stress2/E)^(R-1)));

%%%%%%%%%%%%%%%%%%%%
% condition function
%%%%%%%%%%%%%%%%%%%%
function [c, ceq] = confun(C,minC,maxC)
c = [minC-C; C-maxC];
ceq = [];