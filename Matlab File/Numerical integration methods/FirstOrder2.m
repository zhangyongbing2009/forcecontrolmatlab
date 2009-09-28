function [udot_0, udotdot_0] = FirstOrder2(u_2,u_1,u_0,dt)

udot_0 = (u_0-u_1)/dt;
udotdot_0 = (u_0-2*u_1+u_2)/dt/dt;
    