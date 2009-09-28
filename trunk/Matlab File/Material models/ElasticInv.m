function varargout = ElasticInv(action,MatData,stress,FID)
%ELASTICINV linear-elastic material for force method
% varargout = ElasticInv(action,MatData,stress,FID)
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

% Written: Andreas Schellenberg (andreas.schellenberg@gmx.net)
% Created: 09/09
% Revision: A

% state variables
persistent stressT;

% extract material properties
tag = MatData.tag;      % unique material tag
F   = MatData.F;        % initial elastic flexibility
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
      strainT = F*stressT(:,tag);
      
      varargout = {strainT};
   % ======================================================================
   case 'getTangent'
      varargout = {F};
   % ======================================================================
   case 'getInitialTangent'
      varargout = {F};
   % ======================================================================
   case 'commitState'
      varargout = {0};      
   % ======================================================================
end
