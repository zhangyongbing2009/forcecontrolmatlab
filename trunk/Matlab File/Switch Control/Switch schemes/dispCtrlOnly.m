function [model analysis state] = dispCtrlOnly(MODEL, ANALYSIS, STATE)
% No switching displacement control only
%
%   controlMode: selected mode of control
%               0 = Displacement Control
%               1 = Force Control
%
% Written by Hong Kim
% Created: 10/21/2009
% Last Update: 10/28/09

[model analysis state] = feval(ANALYSIS.schemeDisp, MODEL, ANALYSIS, STATE);
state.controlMode = 0;

