function varargout = Experimental(action,MatData,trialValue)
%EXPERIMENTAL experimental material
% varargout = Experimental(action,MatData,strain)
%
% action  : switch with following possible values
%              'initialize'         initialize internal variables
%              'setTrialStrain'     set the trial strain
%              'setIncrTrialStrain' set the incr trial strain switch ctrl
%              'setTrialStress'     set the trial stress
%              'setIncrTrialStress' set the incr trial stress switch ctrl
%              'getTrialStrain'     get the current strain that was set
%              'getTrialStress'     get the current stress that was set
%              'getStrain'          get the current strain that was calc
%              'getStress'          get the current stress that was calc
%              'getTangentK'        get the current tangent stiffness
%              'getTangentF'        get the current tangent flexibility
%              'getInitialTangentK' get the initial tangent stiffness
%              'getInitialTangentF' sget the initial tangent flexibility
%              'commitState'        commit state of internal variables
%              'disconnect'         disconnect from experimental site
% MatData : data structure with material information
% strain  : trial strain

% Written: Andreas Schellenberg (andreas.schellenberg@gmx.net)
% Created: 03/07
% Revision: A

% socket identifier
persistent socketID;
% state variables
persistent strainT;
persistent stressT;

persistent FIDd;
persistent FIDf;

% extract material properties
tag    = MatData.tag;            % unique material tag
E      = MatData.E;              % initial elastic modulus
ipAddr = MatData.ipAddr;         % ip-address of simAppSiteServer
ipPort = MatData.ipPort;         % ip-port of simAppSiteServer
if isfield(MatData,'id')
   ndf = length(MatData.id);     % number of degrees of freedom
else
   ndf = 1;
end
if isfield(MatData,'dataSize')
   dataSize = MatData.dataSize;  % size of send and receive vectors
else
   dataSize = 256;
end

% initialize send vector
sData = zeros(1,dataSize);

switch action
   % ======================================================================
   case 'initialize'
      socketID(tag) = TCPSocket('openConnection',ipAddr,ipPort);
      if (socketID(tag)<0)
         error('TCPSocket:openConnection',['Unable to setup connection to ',...
            ipAddr,' : ',num2str(ipPort)]);
      end

      % set the data size for the experimental site
      dataSizes = int32([ndf 0 0 0 0, 0 0 0 ndf 0, dataSize]);
      TCPSocket('sendData',socketID(tag),dataSizes,11);
      
      strainT(:,tag) = zeros(ndf,1);
      stressT(:,tag) = zeros(ndf,1);
      
      % send trial response to experimental site
      sData(1) = 3;
      sData(2:ndf+1) = strainT(:,tag)';
      TCPSocket('sendData',socketID(tag),sData,dataSize);

      % get measured resisting forces
      sData(1) = 10;
      TCPSocket('sendData',socketID(tag),sData,dataSize);
      rData = TCPSocket('recvData',socketID(tag),dataSize);
      stressT(:,tag) = rData(1:ndf)';
      
      FIDd = fopen(['ElementDisp',num2str(tag),'.txt'],'w+');
      FIDf = fopen(['ElementForce',num2str(tag),'.txt'],'w+');
      
      varargout = {stressT(:,tag)};
   % ======================================================================
   case 'setTrialStrain'
      strainT(:,tag) = trialValue;
      
      % send trial response to experimental site
      sData(1) = 3;
      sData(2:ndf+1) = strainT(:,tag)';
      TCPSocket('sendData',socketID(tag),sData,dataSize);
      fprintf(FIDd,'%f\n',trialValue);
      
      varargout = {0};
   % ======================================================================
   case 'setIncrTrialStrain'
      strainT(:,tag) = trialValue;
      
      % send trial response to experimental site
      sData(1) = 3;
      sData(2:ndf+1) = strainT(:,tag)';
      TCPSocket('sendData',socketID(tag),sData,dataSize);
      fprintf(FIDd,'%f\n',strainT(:,tag));
      
      varargout = {0};
   % ======================================================================
   case 'setTrialStress'
      if ~isequal(trialValue,stressT(:,tag))
         stressT(:,tag) = trialValue;
      
         % send trial response to experimental site
         sData(1) = 3;
         sData(2:ndf+1) = stressT(:,tag)';
         TCPSocket('sendData',socketID(tag),sData,dataSize);
      end
      fprintf(FIDf,'%f\n',trialValue);
      
      varargout = {0};
   % ======================================================================
   case 'setIncrTrialStress'
      if ~isequal(trialValue,stressT(:,tag))
         stressT(:,tag) = trialValue;
      
         % send trial response to experimental site
         sData(1) = 3;
         sData(2:ndf+1) = stressT(:,tag)';
         TCPSocket('sendData',socketID(tag),sData,dataSize);
      end
      fprintf(FIDf,'%f\n',stressT(:,tag));
      
      varargout = {0};
   % ======================================================================
   case 'getTrialStain'
      varargout = {strainT(:,tag)};
   % ======================================================================
   case 'getTrialStress'
      varargout = {stressT(:,tag)};
   % ======================================================================
   case 'getStress'
      % get measured resisting forces
      sData(1) = 10;
      TCPSocket('sendData',socketID(tag),sData,dataSize);
      rData = TCPSocket('recvData',socketID(tag),dataSize);
      stressT(:,tag) = rData(1:ndf);
      
      fprintf(FIDf,'%f\n',stressT(:,tag));
      varargout = {stressT(:,tag)};
   % ======================================================================
   case 'getStrain'
      % get measured displacements
      sData(1) = 7;
      TCPSocket('sendData',socketID(tag),sData,dataSize);
      rData = TCPSocket('recvData',socketID(tag),dataSize);
      strainT(:,tag) = rData(1:ndf);
      
      fprintf(FIDd,'%f\n',strainT(:,tag));
      varargout = {strainT(:,tag)};
   % ======================================================================
   case 'getTangentK'
      varargout = {E};
   % ======================================================================
   case 'getTangentF'
      varargout = {1/E};
   % ======================================================================
   case 'getInitialTangentK'
      varargout = {E};
   % ======================================================================
   case 'getInitialTangentF'
      varargout = {1/E};
   % ======================================================================
   case 'commitState'
      sData(1) = 5;
      TCPSocket('sendData',socketID(tag),sData,dataSize);
      
      varargout = {0};
   % ======================================================================
   case 'disconnect'
      sData(1) = 99;
      TCPSocket('sendData',socketID(tag),sData,dataSize);
      TCPSocket('closeConnection',socketID(tag));

      varargout = {0};
   % ======================================================================
end
