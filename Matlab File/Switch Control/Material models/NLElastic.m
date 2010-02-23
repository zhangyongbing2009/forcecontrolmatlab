function varargout = NLElastic(action,MatData,trialValue)
%NLELASTIC nonlinear-elastic material
% varargout = NLElastic(action,MatData,strain)
%
% action  : switch with following possible values
%              'initialize'         initialize internal variables
%              'setTrialStrain'     set the trial strain
%              'setIncrTrialStrain' set the incr trial strain
%              'setTrialStress'     set the trial stress
%              'setIncrTrialStress' set the incr trial stress
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

% Written: Hong Kim (hongkim@berkeley.edu)
% Created: 11/09
% Revision: A

% state variables
persistent strainT;
persistent stressT;
% trial history variables
persistent stressC;
persistent strainC;

persistent FIDd;
persistent FIDf;

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
      strainT(:,tag) = zeros(ndf,1);
      stressT(:,tag) = zeros(ndf,1);
      strainC(:,tag) = zeros(ndf,1);
      stressC(:,tag) = zeros(ndf,1);
      
      FIDd = fopen(['ElementDisp',num2str(tag),'.txt'],'w+');
      FIDf = fopen(['ElementForce',num2str(tag),'.txt'],'w+');
      
      varargout = {0.0};
   % ======================================================================
   case 'setTrialStrain'
      strainT(:,tag) = trialValue;
      fprintf(FIDd,'%12.8f\n',strainT(:,tag));
      
      varargout = {0};
   % ======================================================================
   case 'setIncrTrialStrain'
      strainT(:,tag) = strainC(:,tag) + trialValue;
      fprintf(FIDd,'%12.8f\n',strainT(:,tag));
      
      varargout = {0};
   % ======================================================================
   case 'setTrialStress'
      stressT(:,tag) = trialValue;
      fprintf(FIDf,'%12.8f\n',stressT(:,tag));
      
      varargout = {0};
   % ======================================================================
   case 'setIncrTrialStress'
      stressT(:,tag) = stressC(:,tag) + trialValue;
      fprintf(FIDf,'%12.8f\n',stressT(:,tag));
      
      varargout = {0};
   % ======================================================================
   case 'getTrialStrain'
      varargout = {strainT(:,tag)};
   % ======================================================================
   case 'getTrialStress'
      varargout = {stressT(:,tag)};
   % ======================================================================
   case 'getStrain'
      strainT(:,tag) = -sign(stressT(:,tag))*log(1-sign(stressT(:,tag))*(stressT(:,tag)/amp));     
      fprintf(FIDd,'%12.8f\n',strainT(:,tag));
      
      varargout = {strainT};
   % ======================================================================
   case 'getStress'
      stressT(:,tag) = amp*sign(strainT(:,tag))*(1-exp(-sign(strainT(:,tag))*strainT(:,tag)));
      fprintf(FIDf,'%12.8f\n',stressT(:,tag));
      
      varargout = {stressT};
   % ======================================================================
   case 'getTangentK'
      tangentK = amp*exp(-sign(strainT(:,tag))*strainT(:,tag));
 
      varargout = {tangentK};
   % ======================================================================
   case 'getTangentF'
      tangentF = 1/(amp-sign(stressT(:,tag))*(stressT(:,tag)/amp));
      
      varargout = {tangentF};
   % ======================================================================
   case 'getInitialTangentK'
       E = amp; 
       
       varargout = {E};
   % ======================================================================
   case 'getInitialTangentF'
      E = amp;
       
      varargout = {1/E};
   % ======================================================================
   case 'commitState'
      stressC(:,tag) = stressT(:,tag);
      strainC(:,tag) = strainT(:,tag);
      
      varargout = {0};      
   % ======================================================================
end
