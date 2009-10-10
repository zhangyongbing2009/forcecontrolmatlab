function varargout = BiLinearHystereticForce(action,MatData,stress)
%BILINEARHYSTERETICFORCE bilinear-hysteretic material for force method
% varargout = BiLinearHystereticForce(action,MatData,stress)
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

persistent FID;

% extract material properties
tag  = MatData.tag;        % unique material tag
E    = MatData.E;          % initial elastic modulus
fy   = MatData.fy;         % yield stress
b    = MatData.b;          % hardening ratio
F0   = (1-b)*fy;           % stress at zero strain
ey   = fy/E;               % yield strain
B0   = stressC-E*strainC;  % constant used to define the initial stiffness line which pass throgh the current stress
Fy1  = (F0-b*B0)/(1-b); % Upper bound yield stress limit
Fy2  = Fy1-2*fy;            % Lower bound yield stress limit
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
       % calculate the strain
       if stressT(:,tag) > Fy1
           strainT = ey1+(stressT(:,tag)-Fy1)/b/E;
       elseif stressT(:,tag) < Fy2
           strainT = ey2+(stressT(:,tag)-Fy2)/b/E;
       else
           strainT = (stressT(:,tag)-B0)/E;
       end

       varargout = {strainT};
   % ======================================================================
   case 'getTangent'
      if stressT(:,tag) >= Fy1 || stressT(:,tag) <= Fy2
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
      stressC(tag) = stressT(tag);
      strainC(tag) = strainT(tag);
      
      varargout = {0};
   % ======================================================================
end
