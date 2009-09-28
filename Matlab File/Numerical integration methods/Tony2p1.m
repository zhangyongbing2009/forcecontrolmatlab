function [udot_0, udotdot_0] = Tony2p1(udotdot_1,udot_1,u_1,u_0,dt)
a = udot_1;
b = udotdot_1;
c = 3/dt^3*(u_0-u_1-a*dt-b/2*dt^2);

udot_0 = a+b*dt+c*dt^2;
udotdot_0 = b+2*c*dt;