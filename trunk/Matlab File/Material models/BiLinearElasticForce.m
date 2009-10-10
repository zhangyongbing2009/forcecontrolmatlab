function varargout = BiLinearElasticForce(action,MatData,stress)
%BILINEARELASTICFORCE bilinear-elastic material for force method
% varargout = BiLinearElasticForce(action,MatData,stress)
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

persistent FID;

% extract material properties
tag = MatData.tag;      % unique material tag
E   = MatData.E;        % initial elastic modulus
fy  = MatData.fy;       % yield stress
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
      strainT(:,tag) = zeros(ndf,1);
      stressC(tag) = zeros(ndf,1);
      strainC(tag) = zeros(ndf,1);
      
      FID = fopen(['ElementForce',num2str(tag),'.txt'],'w+');      
      varargout = {0.0};
   % ======================================================================
   case 'setTrialStress'
      stressT(:,tag) = stress;
      fprintf(FID,'%f\n',stress);
      
      varargout = {0};
   % ======================================================================
   case 'getStress'
      varargout = {stressT(:,tag)};
   % ======================================================================
   case 'getStrain'
      if abs(stressT(:,tag)) > fy
         strainT = sign(stressT(:,tag))*(fy/E+(abs(stressT(:,tag))-fy)/(E*b));
      else
         strainT = stressT(:,tag)/E;
      end
      
      varargout = {strainT};
   % ======================================================================
   case 'getTangent'
      if abs(stressT(:,tag)) > fy
         tangentT = 1/(b*E);
      else
         tangentT = 1/E;
      end
      
      varargout = {tangentT};
   % ======================================================================
   case 'getInitialTangent'
      varargout = {1/E};
   % ======================================================================
   case 'commitState'
      varargout = {0};      
   % ======================================================================
end
