function varargout = NLElasticForce(action,MatData,stress)
%NLELASTICFORCE bilinear-elastic material for force method
% varargout = NLElasticForce(action,MatData,stress)
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

% extract material properties
tag = MatData.tag;      % unique material tag
amp  = MatData.amp;     % amplitude

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
      
      varargout = {0};
   % ======================================================================
   case 'getStress'
      varargout = {stressT(:,tag)};
   % ======================================================================
   case 'getStrain'
      strainT = -sign(stressT(:,tag))*log(1-sign(stressT(:,tag))*(stressT(:,tag)/amp));
      
      varargout = {strainT};
   % ======================================================================
   case 'getTangent'
      tangentT = 1/(amp-sign(stressT(:,tag))*(stressT(:,tag)/amp));
      
      varargout = {tangentT};
   % ======================================================================
   case 'getInitialTangent'
      E = amp;
       
      varargout = {1/E};
   % ======================================================================
   case 'commitState'
      varargout = {0};      
   % ======================================================================
end
