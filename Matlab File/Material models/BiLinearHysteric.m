function varargout = BiLinearHysteric(action,MatData,strain)
%BILINEARHYSTERIC bilinear-elastic material for displacement method
% varargout = BiLinearHysteric(action,MatData,stress)
%
% action  : switch with following possible values
%              'initialize'         initialize internal variables
%              'setTrialStress'     set the trial stress
%              'getStress'          get the current stress
%              'getStrain'          get the current strain
%              'getTangent'         get the current tangent flexibility
%              'getInitialTangent'  get the initial tangent flexibility
%              'commitState'        commit state of internal variables
% MatData : data structure with material information
% stress  : trial stress

% Written: T.Y. Yang (yangtony2004@gmail.com)
% Created: 09/09
% Revision: A

% state variables
persistent stressT;
persistent strainT;
% trial history variables
persistent stressC;
persistent strainC;

% extract material properties
tag  = MatData.tag;        % unique material tag
E    = MatData.E;          % initial elastic modulus
Fy   = MatData.Fy;         % yield stress
b    = MatData.b;          % hardening ratio
F0   = (1-b)*Fy;           % stress at zero strain
ey   = Fy/E;               % yield strain
B0   = stressC-E*strainC;  % constant used to define the initial stiffness line which pass throgh the current stress
Fy1  = (F0-b*B0)/(1-b); % Upper bound yield stress limit
Fy2  = Fy1-2*Fy;            % Lower bound yield stress limit
ey1  = (F0-B0)/E/(1-b); % Upper bound yield strain limit
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
      
      varargout = {0.0};
   % ======================================================================
   case 'setTrialStrain'
      strainT(:,tag) = strain;      
      varargout = {0};
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
       
      varargout = {stressT};
   % ======================================================================
   case 'getStrain'   
      varargout = {strainT(:,tag)};
   % ======================================================================
   case 'getTangent'
      if strainT(:,tag) > ey1 || strainT(:,tag) < ey2
         tangentT = b*E;
      else
         tangentT = E;
      end      
      varargout = {tangentT};
   % ======================================================================
   case 'getInitialTangent'
      varargout = {E};
   % ======================================================================
   case 'commitState'
      stressC(tag) = stressT(tag);
      strainC(tag) = strainT(tag);
      varargout = {0};      
   % ======================================================================
end
