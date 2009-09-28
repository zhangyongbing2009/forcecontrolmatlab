% Newmark equation
function [udot_0, udotdot_0] = linearAcceleration(udotdot_1,udot_1,u_1,u_0,dt)

udotdot_0 = 6/dt/dt*(u_0-u_1-dt*udot_1)-2*udotdot_1;
udot_0    = udot_1+dt/2*(udotdot_1+udotdot_0);
