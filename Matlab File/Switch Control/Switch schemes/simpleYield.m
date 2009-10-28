function [model analysis state] = simpleYield(MODEL, ANALYSIS, STATE)
% This simple function switches between disp and force algor depending
% on the user defined yield strength of the experimental element.
%
%   controlMode: selected mode of control
%               0 = Displacement Control
%               1 = Force Control
%
% Written by Hong Kim
% Created: 10/21/2009
% Last Update: 10/28/09

threshold = MODEL.MatData(1).fy;
p  = STATE.pr(1,1); 

if (abs(p) <= threshold) 
    [model analysis state] = feval(ANALYSIS.schemeForce, MODEL, ANALYSIS, STATE);
    state.controlMode = 1;
else
    [model analysis state] = feval(ANALYSIS.schemeDisp, MODEL, ANALYSIS, STATE);
    state.controlMode = 0;
end
