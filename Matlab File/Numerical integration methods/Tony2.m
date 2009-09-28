function [udot_0, udotdot_0] = Tony2(udotdot_1,udot_1,u_1,u_0,dt)
A = [1 -dt dt^2; dt -dt^2/2 dt^3/3; 0 1 -2*dt];
B = [udot_1; u_0-u_1; udotdot_1];
X = A\B;
a = X(1);
b = X(2);
c = X(3);

udot_0 = a;
udotdot_0 = b;