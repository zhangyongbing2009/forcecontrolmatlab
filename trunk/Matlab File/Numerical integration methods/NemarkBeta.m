% Newmark equation
function [udot_0, udotdot_0] = NemarkBeta(udotdot_1,udot_1,u_1,u_0,beta,dt)

udotdot_0 = (u_0-u_1-dt*udot_1)/dt/dt/beta-(1/2-beta)*udotdot_1/beta;
udot_0    = udot_1+dt/2*(udotdot_1+udotdot_0);
