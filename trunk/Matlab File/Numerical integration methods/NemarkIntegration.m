% Newmark equation
function [udot_0, udotdot_0] = NemarkIntegration(udotdot_1,udot_1,u_1,u_0,gamma,beta,dt)
c2 = gamma/beta/dt;
c3 = 1/beta/dt/dt;
a1 = 1-gamma/beta;
a2 = dt*(1-1/2*gamma/beta);
a3 = -1/beta/dt;
a4 = 1-1/2/beta;

udot_0 = c2*(u_0-u_1)+a1*udot_1+a2*udotdot_1;
udotdot_0 = c3*(u_0-u_1)+a3*udot_1+a4*udotdot_1;