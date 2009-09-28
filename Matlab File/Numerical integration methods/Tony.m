function [udot_0, udotdot_0] = Tony(udotdot_1,udot_2,udot_1,u_1,u_0,dt)
a = udot_1;
b = udotdot_1;
A = [dt^2 -dt^3; dt^3/3 dt^4/4];
B = [udot_2-a+b*dt;  u_0-u_1-a*dt-b/2*dt^2];
X = A\B;
c = X(1);
d = X(2);

udot_0 = a+b*dt+c*dt^2+d*dt^3;
udotdot_0 = b+2*c*dt+3*d*dt^2;