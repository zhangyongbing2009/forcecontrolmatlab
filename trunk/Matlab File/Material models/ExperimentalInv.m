function varargout = ExperimentalInv(action,MatData,stress)
%EXPERIMENTALINV experimental material for force method
% varargout = ExperimentalInv(action,MatData,stress)
%
% action  : switch with following possible values
%              'initialize'         initialize internal variables
%              'setTrialStress'     set the trial stress
%              'getStress'          get the current stress
%              'getStrain'          get the current strain
%              'getTangent'         get the current tangent flexibility
%              'getInitialTangent'  get the initial tangent flexibility
%              'commitState'        commit state of internal variables
%              'disconnect'         disconnect from experimental site
% MatData : data structure with material information
% stress  : trial stress

% Written: Andreas Schellenberg (andreas.schellenberg@gmx.net)
% Created: 03/07
% Revision: A

% socket identifier
persistent socketID;
% state variables
persistent stressT;

% extract material properties
tag    = MatData.tag;            % unique material tag
F      = MatData.F;              % initial elastic flexibility
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
      dataSizes = int32([0 0 0 ndf 0, ndf 0 0 0 0, dataSize]);
      TCPSocket('sendData',socketID(tag),dataSizes,11);
      
      stressT(:,tag) = zeros(ndf,1);
      
      % send trial response to experimental site
      sData(1) = 3;
      sData(2:ndf+1) = stressT(:,tag)';
      TCPSocket('sendData',socketID(tag),sData,dataSize);

      % get measured displacements
      sData(1) = 7;
      TCPSocket('sendData',socketID(tag),sData,dataSize);
      rData = TCPSocket('recvData',socketID(tag),dataSize);
      strainT = rData(1:ndf)';
      
      varargout = {strainT};
   % ======================================================================
   case 'setTrialStress'
      stressT(:,tag) = stress;
      
      % send trial response to experimental site
      sData(1) = 3;
      sData(2:ndf+1) = stressT(:,tag)';
      TCPSocket('sendData',socketID(tag),sData,dataSize);
      
      varargout = {0};
   % ======================================================================
   case 'getStress'
      varargout = {stressT(:,tag)};
   % ======================================================================
   case 'getStrain'
      % get measured displacements
      sData(1) = 7;
      TCPSocket('sendData',socketID(tag),sData,dataSize);
      rData = TCPSocket('recvData',socketID(tag),dataSize);
      strainT = rData(1:ndf);
      
      varargout = {strainT};
   % ======================================================================
   case 'getTangent'
      varargout = {F};
   % ======================================================================
   case 'getInitialTangent'
      varargout = {F};
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
