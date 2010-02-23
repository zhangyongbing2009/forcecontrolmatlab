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

fy = MODEL.MatData(1).fy;
Rd = .95*fy;
Rf = .9*fy;
p  = STATE.pr(1,1); 
controlMode = STATE.controlMode;
controlModeNew = controlMode;

if (controlMode == 0)
   if (abs(p) < Rf) 
      controlModeNew = 1;
   end
elseif (controlMode == 1)
   if (abs(p) > Rd) 
      controlModeNew = 0;
   end
end
   
if (controlModeNew == 0)
   [model analysis state] = feval(ANALYSIS.schemeDisp, MODEL, ANALYSIS, STATE);
elseif (controlModeNew ==1)
   [model analysis state] = feval(ANALYSIS.schemeForce, MODEL, ANALYSIS, STATE);
end

state.controlMode = controlModeNew;
state.k = 1;
