function varargout = BiLinearElastic(action,MatData,trialValue)
%BILINEARELASTIC bilinear-elastic material
% varargout = BiLinearElastic(action,MatData,trialValue)
%
% action  : switch with following possible values
%              'initialize'         initialize internal variables
%              'setTrialStrain'     set the trial strain
%              'setIncrTrialStrain' set the trial incr strain
%              'setTrialStress'     set the trial stress
%              'setIncrTrialStress' set the trial incr stress
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

% Written: Tony Yang (yangtony2004@gmail.com)
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
      stressC(:,tag) = zeros(ndf,1);
      strainC(:,tag) = zeros(ndf,1);
      
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
      if abs(stressT(:,tag)) > fy
         strainT = sign(stressT(:,tag))*(fy/E+(abs(stressT(:,tag))-fy)/(E*b));
      else
         strainT = stressT(:,tag)/E;
      end
      
      fprintf(FIDd,'%12.8f\n',strainT);
      varargout = {strainT(:,tag)};
   % ======================================================================
   case 'getIncrStrain'
      if abs(stressT(:,tag)) > fy
         strainT = sign(stressT(:,tag))*(fy/E+(abs(stressT(:,tag))-fy)/(E*b));
      else
         strainT = stressT(:,tag)/E;
      end
      
      fprintf(FIDd,'%12.8f\n',strainT);
      varargout = {strainT(:,tag)-strainC(:,tag)};
   % ======================================================================
   case 'getStress'
      if abs(strainT(:,tag)) > fy/E
         stressT = sign(strainT(:,tag))*(fy+(abs(strainT(:,tag))-fy/E)*(E*b));
      else
         stressT = strainT(:,tag)*E;
      end
      
      fprintf(FIDf,'%12.8f\n',stressT);
      varargout = {stressT(:,tag)};
   % ======================================================================
   case 'getIncrStress'
      if abs(strainT(:,tag)) > fy/E
         stressT = sign(strainT(:,tag))*(fy+(abs(strainT(:,tag))-fy/E)*(E*b));
      else
         stressT = strainT(:,tag)*E;
      end
      
      fprintf(FIDf,'%12.8f\n',stressT);
      varargout = {stressT(:,tag)-stressC(:,tag)};
   % ======================================================================
   case 'getTangentK'
      if abs(strainT(:,tag)) >= fy/E
         tangentT = b*E;
      else
         tangentT = E;
      end
      
      varargout = {tangentT};
   % ======================================================================
   case 'getTangentF'
      if abs(stressT(:,tag)) >= fy
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
      stressC(:,tag) = stressT(:,tag);
      strainC(:,tag) = strainT(:,tag);
      
      varargout = {0};      
   % ======================================================================
end
