function varargout = HardeningForce(action,MatData,stress)
%HARDENINGFORCE hardening material (one-dimensional rate-independent plasticity
% model with combined isotropic and kinematic hardening) for force method
% varargout = HardeningForce(action,MatData,strain)
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
% Created: 03/07
% Revision: A

% committed history variables
persistent stressPlasticC;
persistent strainBackC; 
% trial history variables
persistent stressPlasticT;
persistent strainBackT;
% state variables
persistent stressT;
persistent strainT;
persistent tangentT;

% extract material properties
tag  = MatData.tag;     % unique material tag
E    = MatData.E;       % initial elastic modulus
Fy   = MatData.Fy;      % yield stress
Hkin = MatData.Hkin;    % kinematic hardening modulus
F    = 1/E;             % initial elastic flexibility
epsy = Fy/E;            % yield strain
Fkin = 1/Hkin;          % kinematic hardening flexibility

switch action
   % ======================================================================
   case 'initialize'
      stressPlasticC(tag) = 0.0;
      strainBackC(tag)    = 0.0;
      
      stressPlasticT(tag) = 0.0;
      strainBackT(tag)    = 0.0;
      
      stressT(tag)  = 0.0;
      strainT(tag)  = 0.0;
      tangentT(tag) = F;
      
      varargout = {stressT(tag)};
   % ======================================================================
   case 'setTrialStress'
      % check for fast return
      if (abs(stressT(tag)-stress) > eps)
         % set trial stress
         stressT(tag) = stress;
         
         % elastic trial strain
         strainT(tag) = F*(stressT(tag) - stressPlasticC(tag));
         
         % compute trial strain relative to committed back strain
         xi = strainT(tag) - strainBackC(tag);
         
         % compute yield criterion
         Y = abs(xi) - epsy;
         
         % elastic step -> no updates required
         if (Y <= -F*eps)
            % set trial tangent
            tangentT(tag) = F;
            
         % plastic step -> return mapping
         else
            % compute consistency parameter
            dGamma = Y/(F+Fkin);
                        
            % find sign of xi
            signXi = sign(xi);

            % bring trial strain back to yield surface
            strainT(tag) = strainT(tag) - dGamma*F*signXi;
            
            % update plastic stress
            stressPlasticT(tag) = stressPlasticC(tag) + dGamma*signXi;
            
            % update back strain
            strainBackT(tag) = strainBackC(tag) + dGamma*Fkin*signXi;
            
            % set trial tangent
            tangentT(tag) = F*Fkin/(F+Fkin);
         end
      end
      
   	varargout = {0};
   % ======================================================================
   case 'getStress'      
      varargout = {stressT(tag)};
   % ======================================================================
   case 'getStrain'      
      varargout = {strainT(tag)};
   % ======================================================================
   case 'getTangent'      
      varargout = {tangentT(tag)};
   % ======================================================================
   case 'getInitialTangent'
      varargout = {F};
   % ======================================================================
   case 'commitState'
      stressPlasticC(tag) = stressPlasticT(tag);
      strainBackC(tag) = strainBackT(tag);
      
      varargout = {0};      
   % ======================================================================
end
