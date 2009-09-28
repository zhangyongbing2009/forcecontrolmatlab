function [udot_0, udotdot_0] = Diff_Hermite(udot_past,u_past,dt,method)

% check the dimension is consistant
if size(u_past,1) > size(u_past,2)
    u_past = u_past';
end
if size(udot_past,1) > size(udot_past,2)
    udot_past = udot_past';
end
switch method
    case 'pchip'
        % form cubic hermite polynomial
        pp = pchip((-(length(u_past)-1):0)*dt,u_past);
    case 'spline'
        % form spline polynomial
        pp = spline((-(length(u_past)-1):0)*dt,u_past);
end
% take derivative
udot_0 = polyval(polyder(pp.coefs(end,:)),0);

udot_past = [udot_past, udot_0];
switch method
    case 'pchip'
        % form cubic hermite polynomial
        pp = pchip((-(length(udot_past)-1):0)*dt,udot_past);
    case 'spline'
        % form spline polynomial
        pp = spline((-(length(udot_past)-1):0)*dt,udot_past);
end
% take derivative
udotdot_0 = polyval(polyder(pp.coefs(end,:)),0);
