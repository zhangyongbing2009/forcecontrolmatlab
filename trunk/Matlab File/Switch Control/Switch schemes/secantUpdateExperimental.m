function [model analysis state] = secantUpdateExperimental(MODEL, ANALYSIS, STATE)
% This function switches between disp and force control using the algor
% used by Elkhoraibi and Mosalam (2007)
%
%   controlMode: selected mode of control
%               0 = Displacement Control
%               1 = Force Control
%
% Written by Hong Kim
% Created: 2/5/2010
% Last Update: 2/5/2010

% Initialize Variables
kd = ANALYSIS.Kd;
kf = ANALYSIS.Kf;
rd = ANALYSIS.Rd;
rf = ANALYSIS.Rf;
controlMode = STATE.controlMode;
controlModeNew = controlMode;
switchFlag = STATE.switchFlag;

r  = STATE.pr(1,1);
rPrev  = STATE.prPrev(1,1); 
u = STATE.u(1,1);
uPrev = STATE.uPrev(1,1);

% Calculate k
k = (r-rPrev)/(u-uPrev);

if (controlMode == 0)
    if (k > kf)
        controlModeNew = 1;
        switchFlag = 1;
    end 
elseif (controlMode == 1)
    if (k < kd)
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
state.k = k;
