function varargout = BiLinearElasticInv(action,MatData,stress,FID)
%BILINEARELASTICINV bilinear-elastic material for force method
% varargout = BiLinearElasticInv(action,MatData,stress,FID)
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

% Written: Tony Yang (yangtony2004@gmail.com)
% Created: 09/09
% Revision: A

% state variables
persistent stressT;

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
      stressT(:,tag) = zeros(ndf,1);
      
      varargout = {0.0};
   % ======================================================================
   case 'setTrialStress'
      stressT(:,tag) = stress;
      if (nargin==4) && ~isempty(FID)
         fprintf(FID,'%f\n',stressT(:,tag));
      end
      
      varargout = {0};
   % ======================================================================
   case 'getStress'
      varargout = {stressT(:,tag)};
   % ======================================================================
   case 'getStrain'
      if abs(stressT(:,tag)) > Fy
         strainT = sign(stressT(:,tag))*(Fy/E+(abs(stressT(:,tag))-Fy)/(E*b));
      else
         strainT = stressT(:,tag)/E;
      end
      
      varargout = {strainT};
   % ======================================================================
   case 'getTangent'
      if abs(stressT(:,tag)) > Fy
         F = 1/(b*E);
      else
         F = 1/E;
      end
      
      varargout = {F};
   % ======================================================================
   case 'getInitialTangent'
      varargout = {1/E};
   % ======================================================================
   case 'commitState'
      varargout = {0};      
   % ======================================================================
end
