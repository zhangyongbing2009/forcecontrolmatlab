function varargout = BiLinearHysteretic(action,MatData,trialValue)
%BILINEARHYSTERETIC bilinear-hysteretic material for displacement method
% varargout = BiLinearHysteretic(action,MatData,trialValue)
%
% action  : switch with following possible values
%              'initialize'         initialize internal variables
%              'setTrialStrain'     set the trial strain
%              'setTrialStress'     set the trial stress
%              'getTrialStrain'     get the current strain that was set
%              'getTrialStress'     get the current stress that was set
%              'getStrain'          get the current strain that was calc
%              'getStress'          get the current stress that was calc
%              'getTangentK'        get the current tangent stiffness
%              'getTangentF'        get the current tangent flexibility
%              'getInitialTangentK' get the initial tangent stiffness
%              'getInitialTangentF' get the initial tangent flexibility
%              'commitState'        commit state of internal variables
% MatData : data structure with material information
% trialValue  : trial strain or stress depending on the action


% Written: T.Y. Yang (yangtony2004@gmail.com)
% Modified: Hong Kim (hongkim@berkeley.edu)
% Created: 10/09
% Revision: B

% state variables
persistent stressT;
persistent strainT;
% trial history variables
persistent stressC;
persistent strainC;

persistent FIDd;
persistent FIDf;

% extract material properties
tag  = MatData.tag;        % unique material tag
E    = MatData.E;          % initial elastic modulus
fy   = MatData.fy;         % yield stress
b    = MatData.b;          % hardening ratio
F0   = (1-b)*fy;           % stress at zero strain
ey   = fy/E;               % yield strain
B0   = stressC-E*strainC;  % constant used to define the initial stiffness line which pass throgh the current stress
Fy1  = (F0-b*B0)/(1-b);    % Upper bound yield stress limit
Fy2  = Fy1-2*fy;           % Lower bound yield stress limit
ey1  = (F0-B0)/E/(1-b);    % Upper bound yield strain limit
ey2  = ey1-2*ey;           % Lower bound yield strain limit

if isfield(MatData,'id')
   ndf = length(MatData.id);     % number of degrees of freedom
else
   ndf = 1;
end

switch action
   % ======================================================================
   case 'initialize'
      stressT(:,tag) = zeros(ndf,1);
      strainT(:,tag) = zeros(ndf,1);
      stressC(tag) = zeros(ndf,1);
      strainC(tag) = zeros(ndf,1);
      
      FIDd = fopen(['ElementDisp',num2str(tag),'.txt'],'w+');
      FIDf = fopen(['ElementForce',num2str(tag),'.txt'],'w+');
      
      varargout = {0.0};
   % ======================================================================
   case 'setTrialStrain'
      strainT(:,tag) = trialValue;
      fprintf(FIDd,'%f\n',trialValue);
      
      varargout = {0};
   % ======================================================================
   case 'setTrialStress'
      stressT(:,tag) = trialValue;
      fprintf(FIDf,'%f\n',trialValue);
      
      varargout = {0};
   % ======================================================================
   case 'getTrialStrain'   
      varargout = {strainT(:,tag)};
   % ======================================================================
   case 'getTrialStress'
      varargout = {stressT(:,tag)};
   % ======================================================================
   case 'getStrain'
       % calculate the strain
       if stressT(:,tag) > Fy1
           strainT = ey1+(stressT(:,tag)-Fy1)/b/E;
       elseif stressT(:,tag) < Fy2
           strainT = ey2+(stressT(:,tag)-Fy2)/b/E;
       else
           strainT = (stressT(:,tag)-B0)/E;
       end
       
       fprintf(FIDd,'%f\n',strainT);
       varargout = {strainT};
   % ======================================================================
   case 'getStress'
       % calculate the stress
       if strainT(:,tag) > ey1
           stressT = Fy1+(strainT(:,tag)-ey1)*b*E;
       elseif strainT(:,tag) < ey2
           stressT = Fy2+(strainT(:,tag)-ey2)*b*E;
       else
           stressT = strainT(:,tag)*E+B0;
       end
       
       fprintf(FIDf,'%f\n',stressT);
       varargout = {stressT};
   % ======================================================================
   case 'getTangentK'
      if strainT(:,tag) >= ey1 || strainT(:,tag) <= ey2
         tangentT = b*E;
      else
         tangentT = E;
      end
      
      varargout = {tangentT};
   % ======================================================================
   case 'getTangentF'
      if stressT(:,tag) >= Fy1 || stressT(:,tag) <= Fy2
         tangentT = 1/(b*E);
      else
         tangentT = 1/E;
      end
      
      varargout = {tangentT};
   % ======================================================================
   case 'getInitialTangentK'
      varargout = {E};
   % ======================================================================
   case 'getInitialTangentF'
      varargout = {1/E};
   % ======================================================================
   case 'commitState'
      stressC(tag) = stressT(tag);
      strainC(tag) = strainT(tag);
      
      varargout = {0};      
   % ======================================================================
end
