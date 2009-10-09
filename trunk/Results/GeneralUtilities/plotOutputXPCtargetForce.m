function data = plotOutputXPCtargetForce(fileName,actID,iDelay,tStart)
%PLOTOUTPUTXPCTARGETFORCE to plot the xPC-Target output from a hybrid simulation
% data = plotOutputXPCtargetForce(fileName,actID,iDelay,tStart)
%
% data     : data structure
% fileName : .mat file to be loaded
% actID    : id of actuator to plot (optional)
% tStart   : time to start plotting (optional)
%
% Written: Andreas Schellenberg (andreas.schellenberg@gmx.net)
% Created: 04/05

if (nargin<2)
   actID = 1;
   iDelay = 0;
   id = 1024;
elseif (nargin<3)
   iDelay = 0;
   id = 1024;
elseif (nargin<4)
   id = 1024;   
else
   id = tStart*1024;
end

% load the file and extract data
data = [];
%iDelay = 93;
%iDelay = 118;
load(fileName);
targFrc = data(1).values;
commFrc = data(2).values;
measDsp = data(3).values;
state   = data(4).values;
counter = data(5).values;
flag    = data(6).values;
measFrc = data(7).values;

% get screen size
SS = get(0,'screensize');

%==========================================================================
% command forces
CreateWindow('cen',0.80*SS(4)/3*4,0.80*SS(4));
try
   plot(commFrc(id:end,end),commFrc(id:end,actID),'-');
end
grid('on');
xlabel('Time [sec]');
ylabel('Command Force [kip]');
title('Command Force from xPC-Target');
%==========================================================================
% target, command and measured forces
CreateWindow('cen',0.80*SS(4)/3*4,0.80*SS(4));
try
   plot(targFrc(id:end,end),targFrc(id:end,actID),'-b');
end
hold('on');
try
   plot(commFrc(id:end,end),commFrc(id:end,actID),'-r');
end
try
   plot(measFrc(id:end,end),measFrc(id:end,actID),'-g');
end
% try
%    plot(state(id:end,end),state(id:end,1),'-k');
% end
grid('on');
xlabel('Time [sec]');
ylabel('Force [kip]');
title('Forces from xPC-Target');
legend('targFrc','commFrc','measFrc');
%==========================================================================
% error between measured and target forces
CreateWindow('cen',0.80*SS(4)/3*4,0.80*SS(4));
try
   if (iDelay==0)
      N = max(counter(:,1));
      targID = find(counter(:,1) == N);
   else      
      targID = find(counter(:,1) >= 1+iDelay & counter(:,1) < 2+iDelay);
   end
   tID = find(targID >= id);
   targID = targID(tID);
   error = measFrc(targID,actID) - targFrc(targID,actID);
   plot(measFrc(targID,end),targFrc(targID,actID),'-b');
   hold('on');
   plot(measFrc(targID,end),measFrc(targID,actID),'-r');
   plot(measFrc(targID,end),error,'-g');
end
grid('on');
xlabel('Time [sec]');
ylabel('Forces [kip]');
title('Force Errors from xPC-Target');
legend('targFrc','measFrc','error');
%==========================================================================
% fft of error
CreateWindow('cen',0.80*SS(4)/3*4,0.80*SS(4));
try
   N = max(counter(:,1));
   if (iDelay==0)
      targID = find(counter(:,1) == N);
   else      
      targID = find(counter(:,1) >= 1+iDelay & counter(:,1) < 2+iDelay);
   end
   tID = find(targID >= id);
   targID = targID(tID);
   error = measFrc(targID,actID) - targFrc(targID,actID);
   dt = N/1024;
   getFFT(error,dt,'Error between Measured and Target Force');
end
%==========================================================================
% measured vs. target forces
CreateWindow('cen',0.80*SS(4)/3*4,0.80*SS(4));
try
   if (iDelay==0)
      N = max(counter(:,1));
      targID = find(counter(:,1) == N);
   else      
      targID = find(counter(:,1) >= 1+iDelay & counter(:,1) < 2+iDelay);
   end
   tID = find(targID >= id);
   targID = targID(tID);
   plot(targFrc(targID,actID),measFrc(targID,actID),'-b');
   hold('on');
end
grid('on');
xlabel('Target Force [ikip]');
ylabel('Measured Force [kip]');
title('Forces from xPC-Target');
%==========================================================================
% measured displacement
CreateWindow('cen',0.80*SS(4)/3*4,0.80*SS(4));
try
   plot(measDsp(id:end,end),measDsp(id:end,actID),'-');
end
grid('on');
xlabel('Time [sec]');
ylabel('Measured Displacement [in.]');
title('Measured Displacement from xPC-Target');
%==========================================================================
% fft of measured displacement
CreateWindow('cen',0.80*SS(4)/3*4,0.80*SS(4));
try
   dt = 1/1024;
   getFFT(measDsp(id:end,actID),dt,'Measured Displacement');
   set(gca,'YScale','log');
end
%==========================================================================
% hysteresis loop
CreateWindow('cen',0.80*SS(4)/3*4,0.80*SS(4));
try
   if (iDelay==0)
      N = max(counter(:,1));
      targID = find(counter(:,1) == N);
   else      
      targID = find(counter(:,1) >= 1+iDelay & counter(:,1) < 2+iDelay);
   end
   tID = find(targID >= id);
   targID = targID(tID);
   plot(measDsp(targID,actID),measFrc(targID,actID),'-');
end
grid('on');
xlabel('Measured Displacement [in.]');
ylabel('Measured Force [kip]');
title('Hysteresis Loop from xPC-Target');
%==========================================================================
% state of predictor-corrector
CreateWindow('cen',0.80*SS(4)/3*4,0.80*SS(4));
try
   plot(state(id:end,end),state(id:end,1),'-b');
end
grid('on');
xlabel('Time [sec]');
ylabel('State [-]');
title('State from xPC-Target');
%==========================================================================
% counter of predictor-corrector
CreateWindow('cen',0.80*SS(4)/3*4,0.80*SS(4));
try
   plot(counter(id:end,end),counter(id:end,1),'-b');
end
grid('on');
xlabel('Time [sec]');
ylabel('Counter [-]');
title('Counter from xPC-Target');
%==========================================================================
% update and target flags of predictor-corrector
CreateWindow('cen',0.80*SS(4)/3*4,0.80*SS(4));
try
   plot(flag(id:end,end),flag(id:end,1),'-b');
end
hold('on');
try
   plot(flag(id:end,end),flag(id:end,2),'-r');
end
try
   plot(flag(id:end,end),flag(id:end,3),'-g');
end
grid('on');
xlabel('Time [sec]');
ylabel('Flag [-]');
title('Flags from xPC-Target');
legend('newTarget','switchPC','atTarget');
%==========================================================================
