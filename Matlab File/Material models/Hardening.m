function varargout = Hardening(action,MatData,strain)
%HARDENING hardening material (one-dimensional rate-independent plasticity
% model with combined isotropic and kinematic hardening)
% varargout = Hardening(action,MatData,strain)
%
% action  : switch with following possible values
%              'initialize'         initialize internal variables
%              'setTrialStrain'     set the trial strain
%              'getStrain'          get the current strain
%              'getStress'          get the current stress
%              'getTangent'         get the current tangent modulus
%              'getInitialTangent'  get the initial tangent modulus
%              'commitState'        commit state of internal variables
% MatData : data structure with material information
% strain  : trial strain

% Written: Andreas Schellenberg (andreas.schellenberg@gmx.net)
% Created: 03/07
% Revision: A

% committed history variables
persistent strainPlasticC;
persistent stressBackC; 
persistent hardeningC;
% trial history variables
persistent strainPlasticT;
persistent stressBackT;
persistent hardeningT;
% state variables
persistent strainT;
persistent stressT;
persistent tangentT;

% extract material properties
tag  = MatData.tag;     % unique material tag
E    = MatData.E;       % initial elastic modulus
fy   = MatData.fy;      % yield stress
Hkin = MatData.Hkin;    % kinematic hardening modulus
Kiso = MatData.Kiso;    % isotropic hardening modulus

switch action
   % ======================================================================
   case 'initialize'
      strainPlasticC(tag) = 0.0;
      stressBackC(tag)    = 0.0;
      hardeningC(tag)     = 0.0;
      
      strainPlasticT(tag) = 0.0;
      stressBackT(tag)    = 0.0;
      hardeningT(tag)     = 0.0;
      
      strainT(tag)  = 0.0;
      stressT(tag)  = 0.0;
      tangentT(tag) = E;
      
      varargout = {stressT(tag)};
   % ======================================================================
   case 'setTrialStrain'
      % check for fast return
      if (abs(strainT(tag)-strain) > eps)
         % set trial strain
         strainT(tag) = strain;
         
         % elastic trial stress
         stressT(tag) = E*(strainT(tag) - strainPlasticC(tag));
         
         % compute trial stress relative to committed back stress
         xi = stressT(tag) - stressBackC(tag);
         
         % compute yield criterion
         Y = abs(xi) - (fy+Kiso*hardeningC(tag));
         
         % elastic step -> no updates required
         if (Y <= -E*eps)
            % set trial tangent
            tangentT(tag) = E;
            
         % plastic step -> return mapping
         else
            % compute consistency parameter
            dGamma = Y/(E+Hkin+Kiso);
                        
            % find sign of xi
            signXi = sign(xi);

            % bring trial stress back to yield surface
            stressT(tag) = stressT(tag) - dGamma*E*signXi;
            
            % update plastic strain
            strainPlasticT(tag) = strainPlasticC(tag) + dGamma*signXi;
            
            % update back stress
            stressBackT(tag) = stressBackC(tag) + dGamma*Hkin*signXi;
            
            % update internal hardening variable
            hardeningT(tag) = hardeningC(tag) + dGamma;
            
            % set trial tangent
            tangentT(tag) = E*(Hkin+Kiso)/(E+Hkin+Kiso);
         end
      end
      
   	varargout = {0};
   % ======================================================================
   case 'getStrain'      
      varargout = {strainT(tag)};
   % ======================================================================
   case 'getStress'      
      varargout = {stressT(tag)};
   % ======================================================================
   case 'getTangent'      
      varargout = {tangentT(tag)};
   % ======================================================================
   case 'getInitialTangent'
      varargout = {E};
   % ======================================================================
   case 'commitState'
      strainPlasticC(tag) = strainPlasticT(tag);
      stressBackC(tag) = stressBackT(tag);
      hardeningC(tag) = hardeningT(tag);
      
      varargout = {0};      
   % ======================================================================
end
