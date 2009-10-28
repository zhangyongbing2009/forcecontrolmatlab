function [model analysis state] = Tangent(MODEL, ANALYSIS, STATE)
% This function determines when to switch between disp and force control
%
% Written by Hong Kim
% Created: 10/21/2009
% Last Update: 10/21/09

[model analysis state] = feval(ANALYSIS.scheme, MODEL, ANALYSIS, STATE);