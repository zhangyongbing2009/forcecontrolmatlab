function varargout = BiLinearElastic(action,MatData,strain)
%BILINEARELASTIC bilinear-elastic material
% varargout = BiLinearElastic(action,MatData,strain)
%
% action  : switch with following possible values
%              'initialize'         initialize internal variables
%              'setTrialStrain'     set the trial strain
%              'getStrain'          get the current strain
%              'getStress'          get the current stress
%              'getTangent'         get the current tangent modulus
%              'getInitialTangent'  get the initial tangent modulus
%              'commitState'        commit state of internal variables
% MatData : data structure with material information
% strain  : trial strain

% Written: Tony Yang (yangtony2004@gmail.com)
% Created: 09/09
% Revision: A

% state variables
persistent strainT;

% extract material properties
tag = MatData.tag;      % unique material tag
E   = MatData.E;        % initial elastic modulus
Fy  = MatData.Fy;       % yield stress
b   = MatData.b;        % hardening ratio

if isfield(MatData,'id')
   ndf = length(MatData.id);     % number of degrees of freedom
else
   ndf = 1;
end

switch action
   % ======================================================================
   case 'initialize'
      strainT(:,tag) = zeros(ndf,1);
      
      varargout = {0.0};
   % ======================================================================
   case 'setTrialStrain'
      strainT(:,tag) = strain;
      
      varargout = {0};
   % ======================================================================
   case 'getStrain'
      varargout = {strainT(:,tag)};
   % ======================================================================
   case 'getStress'
      if abs(strainT(:,tag)) > Fy/E
         stressT = sign(strainT(:,tag))*(Fy+(abs(strainT(:,tag))-Fy/E)*(E*b));
      else
         stressT = strainT(:,tag)*E;
      end
      
      varargout = {stressT};
   % ======================================================================
   case 'getTangent'
      if abs(strainT(:,tag)) > Fy/E
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
      varargout = {0};      
   % ======================================================================
end
