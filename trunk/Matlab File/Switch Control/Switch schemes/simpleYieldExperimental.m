function [model analysis state] = simpleYieldExperimental(MODEL, ANALYSIS, STATE)
% This simple function switches between disp and force algor depending
% on the user defined yield strength of the experimental element.  This
% only works with Elastic Material for now.
%
%   controlMode: selected mode of control
%               0 = Displacement Control
%               1 = Force Control
%
% Written by Hong Kim
% Created: 2/3/2010
% Last Update: 2/3/2010

fy = MODEL.MatData(1).fy;
Rd = ANALYSIS.Rd*fy;
Rf = ANALYSIS.Rf*fy;
p  = STATE.pr(1,1); 
controlMode = STATE.controlMode;
controlModeNew = controlMode;
switchFlag = STATE.switchFlag;

% switching for a softening system
% if (controlMode == 0)
%    if (abs(p) < Rf) 
%       controlModeNew = 1;
%       switchFlag = 1;
%    end
% elseif (controlMode == 1)
%    if (abs(p) > Rd) 
%       controlModeNew = 0;
%       switchFlag = 1;
%    end
% end
   
% switching for a hardening system
if (controlMode == 0)
   if (abs(p) > Rf) 
      controlModeNew = 1;
      switchFlag = 1;
   end
elseif (controlMode == 1)
   if (abs(p) < Rd) 
      controlModeNew = 0;
      switchFlag = 1;
   end
end

% Set the offset so that STS controller can handle the switch only works
% nonhysteric material
if (switchFlag == 1)
   STATE.offsetu = STATE.u;
   STATE.offsetpr = STATE.pr;
   for j=1:MODEL.numElem
      if (isequal(MODEL.Element{j}, 'Experimental') ~= 1)
         feval(MODEL.Element{j},'setTrialStrain',MODEL.MatData(j),STATE.offsetu(j,:));
         feval(MODEL.Element{j},'setTrialStress',MODEL.MatData(j),STATE.offsetpr(j,:));
         feval(MODEL.Element{j},'commitState',MODEL.MatData(j));
      end
   end
end

if (controlModeNew == 0)
   [model analysis state] = feval(ANALYSIS.schemeDisp, MODEL, ANALYSIS, STATE);
elseif (controlModeNew ==1)
   [model analysis state] = feval(ANALYSIS.schemeForce, MODEL, ANALYSIS, STATE);
end

state.offsetu = STATE.offsetu;
state.offsetpr = STATE.offsetpr;
state.controlMode = controlModeNew;
state.switchFlag = switchFlag;
state.k = 1;
