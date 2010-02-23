function [model analysis state] = secantUpdate(MODEL, ANALYSIS, STATE)
% This function switches between disp and force control using the algor
% used by Elkhoraibi and Mosalam (2007)
%
%   controlMode: selected mode of control
%               0 = Displacement Control
%               1 = Force Control
%
% Written by Hong Kim
% Created: 10/21/2009
% Last Update: 10/28/09

% Initialize Variables
kd = ANALYSIS.Kd;
kf = ANALYSIS.Kf;
rd = ANALYSIS.Rd;
rf = ANALYSIS.Rf;
controlMode = STATE.controlMode;
controlModeNew = controlMode;

r  = STATE.pr(1,1);
rPrev  = STATE.prPrev(1,1); 
u = STATE.u(1,1);
uPrev = STATE.uPrev(1,1);

% Calculate k
k = (r-rPrev)/(u-uPrev);

if (controlMode == 0)
    if (k > kf)
        controlModeNew = 1;
    end 
elseif (controlMode == 1)
    if (k < kd)
        controlModeNew = 0;
    end
end

if (controlModeNew == 0)
    [model analysis state] = feval(ANALYSIS.schemeDisp, MODEL, ANALYSIS, STATE);
elseif (controlModeNew ==1)
    [model analysis state] = feval(ANALYSIS.schemeForce, MODEL, ANALYSIS, STATE);
end

state.controlMode = controlModeNew;
state.k = k;
