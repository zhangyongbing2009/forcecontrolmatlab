% constiutive relationship
function [V,DV] = Constutive(Q,n,FID)
    
switch n
    case 1
        Fy = 1.5*1000;
        E = 2.8;
        b = 2;%0.25;
    case 2
        Fy = 10000;
        E = 2;
        b = 0.05;
    case 3
        Fy = 10000;
        E = 5.6;
        b = 0.05;
end
if abs(Q) > Fy
    V = sign(Q)*(Fy/E+(abs(Q)-Fy)/(E*b));
    DV = 1/(b*E);
else
    V = Q/E;
    DV = 1/E;
end

if nargin > 2
    if ~isempty(FID)        
        fprintf(FID,'%f\n',Q);
    end
end
    