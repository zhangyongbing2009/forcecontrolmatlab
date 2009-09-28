function [udot_0, udotdot_0] = SecondOrder(u_3,u_2,u_1,u_0,dt)

udot_0 = 1/2*(3*u_0-4*u_1+u_2)/dt;
udotdot_0 = (2*u_0-5*u_1+4*u_2-u_3)/dt/dt;
    