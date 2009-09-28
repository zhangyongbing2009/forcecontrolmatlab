function varargout = Experimental(action,MatData,strain)
%EXPERIMENTAL experimental material
% varargout = Experimental(action,MatData,strain)
%
% action  : switch with following possible values
%              'initialize'         initialize internal variables
%              'setTrialStrain'     set the trial strain
%              'getStrain'          get the current strain
%              'getStress'          get the current stress
%              'getTangent'         get the current tangent modulus
%              'getInitialTangent'  get the initial tangent modulus
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
      
      % send trial response to experimental site
      sData(1) = 3;
      sData(2:ndf+1) = strainT(:,tag)';
      TCPSocket('sendData',socketID(tag),sData,dataSize);

      % get measured resisting forces
      sData(1) = 10;
      TCPSocket('sendData',socketID(tag),sData,dataSize);
      rData = TCPSocket('recvData',socketID(tag),dataSize);
      stressT = rData(1:ndf)';
      
      varargout = {stressT};
   % ======================================================================
   case 'setTrialStrain'
      strainT(:,tag) = strain;
      
      % send trial response to experimental site
      sData(1) = 3;
      sData(2:ndf+1) = strainT(:,tag)';
      TCPSocket('sendData',socketID(tag),sData,dataSize);
      
      varargout = {0};
   % ======================================================================
   case 'getStain'
      varargout = {strainT(:,tag)};
   % ======================================================================
   case 'getStress'
      % get measured resisting forces
      sData(1) = 10;
      TCPSocket('sendData',socketID(tag),sData,dataSize);
      rData = TCPSocket('recvData',socketID(tag),dataSize);
      stressT = rData(1:ndf);
      
      varargout = {stressT};
   % ======================================================================
   case 'getTangent'
      varargout = {E};
   % ======================================================================
   case 'getInitialTangent'
      varargout = {E};
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
