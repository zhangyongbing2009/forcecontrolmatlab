clear all;
close all;
clc;
% ground motion
path = 'C:\Documents and Settings\Tony\Desktop\Force control\Models\ThreeSpanTruss\';
% get ground-acceleration
ag = load([path,'\elcentro.txt']);
ag = 386.1*ag;

% model parameters
m1 = 0.04;
m2 = 0.02;
c1 = 1.01*m1;
c2 = 1.01*m2;
gamma = 1/2;
beta = 1/4;
dt = 0.02;

data = load([path,'\Node_Dsp.out']);
tN = [0; data(:,1)];
dN = [0 0; data(:,[2,3])];
data = load([path,'\Node_Vel.out']);
vN = [0 0; data(:,[2,3])];
data = load([path,'\Node_Acc.out']);
aN = [0 0; data(:,[2,3])];
% get element data
data = load([path,'\Elmt_Frc.out']);
fE = [0 0 0; data(:,2:4)];
% get element data
dataGF = load([path,'\Elmt_GFrc.out']);
fGE = [0 0 0 0; dataGF(:,[4,6,8,10])];
data = load([path,'\Elmt_Defo.out']);
dE = [0 0 0; data(:,2:4)];

Q1O = fE(:,1); % tension positive
Q2O = fE(:,2); % tension positive
Q3O = fE(:,3); % tension positive
Q1GO = fGE(:,1); % tension positive
Q2GLO = fGE(:,2); % tension positive
Q2GRO = fGE(:,3); % tension positive
Q3GO = fGE(:,4); % tension positive

V1O = dE(:,1); % tension positive
V2O = dE(:,2); % tension positive
V3O = dE(:,3); % tension positive

U1O = dN(:,1);
U2O = dN(:,2);
U1dotO = vN(:,1);
U2dotO = vN(:,2);
U1dotdotO = aN(:,1);
U2dotdotO = aN(:,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check the Nuerical algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
checkFlag = 0;
if checkFlag
    U1Cal = U1O;
    U1dotCal(1:4) = 0;
    U1dotdotCal(1:4) = 0;
    for tt = 5:length(U1Cal)
        %[U1dotCal(tt),U1dotdotCal(tt)] = NemarkIntegration(U1dotdotCal(tt-1),U1dotCal(tt-1),U1Cal(tt-1),U1Cal(tt),gamma,beta,dt);
        %[U1dotCal(tt), U1dotdotCal(tt)] = FirstOrder(U1dotCal(tt-1),U1Cal(tt-1),U1Cal(tt),dt);
        %[U1dotCal(tt), U1dotdotCal(tt)] = FirstOrder2(U1Cal(tt-2),U1Cal(tt-1),U1Cal(tt),dt);
        %[U1dotCal(tt), U1dotdotCal(tt)] = SecondOrder(U1Cal(tt-3),U1Cal(tt-2),U1Cal(tt-1),U1Cal(tt),dt);
        %[U1dotCal(tt), U1dotdotCal(tt)] = ThirdOrder(U1Cal(tt-4),U1Cal(tt-3),U1Cal(tt-2),U1Cal(tt-1),U1Cal(tt),dt);
        %[U1dotCal(tt), U1dotdotCal(tt)] = Diff_Hermite(U1dotCal(tt-4:tt-1),U1Cal(tt-4:tt),dt,'pchip');
        %[U1dotCal(tt), U1dotdotCal(tt)] = Tony(U1dotdotCal(tt-1),U1dotCal(tt-2),U1dotCal(tt-1),U1Cal(tt-1),U1Cal(tt),dt);
        %[U1dotCal(tt), U1dotdotCal(tt)] = Tony2(U1dotdotCal(tt-1),U1dotCal(tt-1),U1Cal(tt-1),U1Cal(tt),dt);
        %[U1dotCal(tt), U1dotdotCal(tt)] = Tony2p1(U1dotdotCal(tt-1),U1dotCal(tt-1),U1Cal(tt-1),U1Cal(tt),dt);
    end
    figure;
    subplot(2,1,1)
    plot(U1dotCal)
    hold on
    plot(U1dotO,'r:')
    ylabel('Udot1')
    grid
    subplot(2,1,2)
    plot(U1dotdotCal)
    hold on
    plot(U1dotdotO,'r:')
    ylabel('Udotdot1')
    grid
end

% equilibrium EQ1
% EQ1O = -m2*ag-m2*aN(:,2)-c2*vN(:,2)+Q3O-Q2O;
% EQ2O = -m1*ag-m1*aN(:,1)-c1*vN(:,1)+Q2O-Q1O;
% % EQ1OG = -m2*ag-m2*aN(:,2)-c2*vN(:,2)-Q3GO-Q2GRO;
% % EQ2OG = -m1*ag-m1*aN(:,1)-c1*vN(:,1)-Q3GO-Q2GRO;
% figure; plot(EQ1O)
% figure; plot(EQ2O)



% figure;
% subplot(3,1,1);
% plot(V1O);
% hold on
% plot(dN(:,1),':r')
% ylabel('V1')
% grid
% subplot(3,1,2);
% plot(V2O);
% hold on
% plot(dN(:,2)-dN(:,1),':r')
% ylabel('V2')
% grid
% subplot(3,1,3);
% plot(V3O);
% hold on
% plot(-dN(:,2),':r')
% ylabel('V3')
% grid
% 
% figure; 
% subplot(3,1,1)
% plot(V1O,Q1O);
% grid
% subplot(3,1,2)
% plot(V2O,Q2O);
% grid
% subplot(3,1,3)
% plot(V3O,Q3O);
% grid

% check the constraints
Q1 = zeros(length(Q3O),1);
Q2 = zeros(length(Q3O),1);
Q3 = zeros(length(Q3O),1);
V1 = zeros(length(Q3O),1);
V2 = zeros(length(Q3O),1);
V3 = zeros(length(Q3O),1);
U1 = zeros(length(Q3O),1);
U2 = zeros(length(Q3O),1);
U1dot = zeros(length(Q3O),1);
U2dot = zeros(length(Q3O),1);
U1dotdot = zeros(length(Q3O),1);
U2dotdot = zeros(length(Q3O),1);

for pp = 5:length(Q3O)
    Q3(pp) = Q3O(pp);
    V3(pp) = Constutive(Q3(pp),3);
    U2(pp) = -V3(pp);
    %[U2dot(pp), U2dotdot(pp)] = NemarkIntegration(U2dotdot(pp-1),U2dot(pp-1),U2(pp-1),U2(pp),gamma,beta,dt);
    %[U2dot(pp), U2dotdot(pp)] = NemarkBeta(U2dotdot(pp-1),U2dot(pp-1),U2(pp-1),U2(pp),beta,dt);
    %[U2dot(pp), U2dotdot(pp)] = linearAcceleration(U2dotdot(pp-1),U2dot(pp-1),U2(pp-1),U2(pp),dt);
    %[U2dot(pp), U2dotdot(pp)] = FirstOrder(U2dot(pp-1),U2(pp-1),U2(pp),dt);
    [U2dot(pp), U2dotdot(pp)] = SecondOrder(U2(pp-3),U2(pp-2),U2(pp-1),U2(pp),dt);
    %[U2dot(pp), U2dotdot(pp)] = ThirdOrder(U2(pp-4),U2(pp-3),U2(pp-2),U2(pp-1),U2(pp),dt);
    %[U2dot(pp), U2dotdot(pp)] = Diff_Hermite(U2dot(pp-4:pp-1),U2(pp-4:pp),dt,'spline');
    
%     clc;
%     U2dot(pp) - vN(pp,2)
%     U2dotdot(pp) - aN(pp,2)
%     
    
    %Q2(pp) = -m2*ag(pp)-m2*U2dotdot(pp)-c2*U2dot(pp)+Q3(pp);
    Q2(pp) = Q2O(pp);
    
    
    V2(pp) = Constutive(Q2(pp),2);
    U1(pp) = U2(pp)-V2(pp);
    %[U1dot(pp), U1dotdot(pp)] = NemarkIntegration(U1dotdot(pp-1),U1dot(pp-1),U1(pp-1),U1(pp),gamma,beta,dt);
    %[U1dot(pp), U1dotdot(pp)] = NemarkBeta(U1dotdot(pp-1),U1dot(pp-1),U1(pp-1),U1(pp),beta,dt);
    %[U1dot(pp), U1dotdot(pp)] = linearAcceleration(U1dotdot(pp-1),U1dot(pp-1),U1(pp-1),U1(pp),dt);
    %[U1dot(pp), U1dotdot(pp)] = FirstOrder(U1dot(pp-1),U1(pp-1),U1(pp),dt);
    [U1dot(pp), U1dotdot(pp)] = SecondOrder(U1(pp-3),U1(pp-2),U1(pp-1),U1(pp),dt);
    %[U1dot(pp), U1dotdot(pp)] = ThirdOrder(U1(pp-4),U1(pp-3),U1(pp-2),U1(pp-1),U1(pp),dt);
    %[U1dot(pp), U1dotdot(pp)] = Diff_Hermite(U1dot(pp-4:pp-1),U1(pp-4:pp),dt,'spline');
    
    Q1(pp) = -m1*ag(pp)-m1*U1dotdot(pp)-c1*U1dot(pp)+Q2(pp);
    V1(pp) = Constutive(Q1(pp),1);
end



% % equilibrium EQs
% EQ1 = -m2*ag-m2*U2dotdot-c2*U2dot+Q3-Q2;
% EQ2 = -m1*ag-m1*U1dotdot-c1*U1dot+Q2-Q1;
% figure; plot(EQ1)
% figure; plot(EQ2)

% plot the figure
figure;
subplot(3,1,1)
plot(Q1(1:pp))
hold on
plot(Q1O(1:pp),':r')
ylabel('Q1')
grid
subplot(3,1,2)
plot(Q2(1:pp))
hold on
plot(Q2O(1:pp),':r')
ylabel('Q2')
grid
subplot(3,1,3)
plot(Q3(1:pp))
hold on
plot(Q3O(1:pp),':r')
ylabel('Q3')
xlabel('Step')
grid

figure;
subplot(3,2,1)
plot(U1(1:pp))
hold on
plot(U1O(1:pp),':r')
ylabel('U1')
grid
subplot(3,2,2)
plot(U2(1:pp))
hold on
ylabel('U2')
plot(U2O(1:pp),':r')
grid
subplot(3,2,3)
plot(U1dot(1:pp))
hold on
plot(U1dotO(1:pp),':r')
ylabel('U1dot')
grid
subplot(3,2,4)
plot(U2dot(1:pp))
hold on
ylabel('U2dot')
plot(U2dotO(1:pp),':r')
grid
subplot(3,2,5)
plot(U1dotdot(1:pp))
hold on
plot(U1dotdotO(1:pp),':r')
ylabel('U1dotdot')
xlabel('Step')
grid
subplot(3,2,6)
plot(U2dotdot(1:pp))
hold on
ylabel('U2dotdot')
xlabel('Step')
plot(U2dotdotO(1:pp),':r')
grid
