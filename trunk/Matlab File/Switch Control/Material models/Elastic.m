function varargout = Elastic(action,MatData,trialValue)
%ELASTIC linear-elastic material
% varargout = Elastic(action,MatData,trialValue)
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
%              'getInitialTangentF' sget the initial tangent flexibility
%              'commitState'        commit state of internal variables
% MatData : data structure with material information
% trialValue  : trial strain or stress depending on the action

% Written: Andreas Schellenberg (andreas.schellenberg@gmx.net)
% Modified: Hong Kim (hongkim@berkeley.edu)
% Created: 10/09
% Revision: B

% state variables
persistent strainT;
persistent stressT;

persistent FIDd;
persistent FIDf;

% extract material properties
tag = MatData.tag;      % unique material tag
E   = MatData.E;        % initial elastic modulus
if isfield(MatData,'id')
   ndf = length(MatData.id);     % number of degrees of freedom
else
   ndf = 1;
end

switch action
   % ======================================================================
   case 'initialize'
      strainT(:,tag) = zeros(ndf,1);
      stressT(:,tag) = zeros(ndf,1);
      
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
      strainT(:,tag) = 1/E*stressT(:,tag);
      fprintf(FIDd,'%f\n',strainT(:,tag));
      varargout = {strainT(:,tag),stressT(:,tag)};
   % ======================================================================
   case 'getStress'
      stressT(:,tag) = E*strainT(:,tag);
      fprintf(FIDf,'%f\n',stressT(:,tag));
      varargout = {stressT(:,tag)};
   % ======================================================================
   case 'getTangentK'
      varargout = {E};
   % ======================================================================
   case 'getInitialTangentK'
      varargout = {E};
   % ======================================================================
   case 'getTangentF'
      varargout = {1/E};
   % ======================================================================
   case 'getInitialTangentF'
      varargout = {1/E};
   % ======================================================================
   case 'commitState'
      varargout = {0};      
   % ======================================================================
end
