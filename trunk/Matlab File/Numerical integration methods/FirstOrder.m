function [udot_0, udotdot_0] = FirstOrder(udot_1,u_1,u_0,dt)

udot_0 = (u_0-u_1)/dt;
udotdot_0 = (udot_0-udot_1)/dt;
    