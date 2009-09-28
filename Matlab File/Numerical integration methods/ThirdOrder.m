function [udot_0, udotdot_0] = ThirdOrder(u_4,u_3,u_2,u_1,u_0,dt)

udot_0 = 1/6*(11*u_0-18*u_1+9*u_2-2*u_3)/dt;
udotdot_0 = 1/12*(35*u_0-104*u_1+114*u_2-56*u_3+11*u_4)/dt/dt;
    