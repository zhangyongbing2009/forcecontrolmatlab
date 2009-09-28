function strain = RambergOsgood_fun(P,K,A,R)
% This is a Matlab implemenation of the Ramberg-Osgood steel material model
% that takes stress/force and returns strain/displacement.  (Matzen and
% McNiven, EERC Report 76-20, 1976, p 97-106)

% IBS number key
% IBS =  1: the ascending skeletal curve
% IBS =  2: the descending skeletal curve
% IBS =  3: the minimum bounding curve
% IBS =  4: the interior ascending branch curve originating from IBS = 2
% IBS =  5: the maximum bounding curve
% IBS =  6: an interior descending branch curve originating from IBS = 1
% IBS =  7: an interior ascending branch curve originating from IBS = 6
% IBS =  8: an interior descending branch curve originating from IBS = 4
% IBS =  9: an interior ascending branch curve originating from IBS = 11
% IBS = 10: an interior descending branch curve originating from IBS = 12
% IBS = 11: an interior descending branch curve originating from IBS = 3
% IBS = 12: an interior ascending branch curve originating from IBS = 5
% IBS = 13: an arbitrary interior ascending branch curve that may originate 
%           from IBS = 8, 10 or 14
% IBS = 14: an arbitrary interior descending branch curve that may 
%           originate from IBS = 7, 9 or 13

% x-Coordinates of zero-force crossing of a curve that is a candidate for a
% bounding curve
%
% strain = stress/E + K*(stress/E)^n
% strain = stress/K*(1 + A*(stress/K)^(R-1)) 

xZ = 0;

% Coordinates of the origin of the min bounding curve
xMINL = 0;  PMINL = 0;

% Coordinates of the terminus of the min bounding curve
xMINU = 0;  PMINU = 0;

% x-Coordinates of the zero-force crossing of the min bounding curve
xMINZ = 0;

% Coordinates of the origin of the max bounding curve
xMAXL = 0;  PMAXL = 0;

% Coordinates of the terminus of the max bounding curve
xMAXU = 0;  PMAXU = 0;

% x-Coordinates of the zero-force crossing of the max bounding curve
xMAXZ = 0;

% Coordinates of the most recent reversal from a bounding or skeletal curve
xBC = 0;    PBC = 0;

% x-Coordinates of the IBS=3 curve corresponding to P, of the current point
x3 = 0;

% x-Coordinates of the IBS=5 curve corresponding to P, of the current point
x5 = 0;

% x-Coordinates of the skeletal curve corresponding to P, of the current pt
xSKEL = 0;

% set parameters K, A, R

n = length(P);  
x = zeros(1,n);
% K = 20;  A = 0.5;  R = 10;

% determine the skeletal curve IBS number
if (P(2) > 0)
    IBS = 1;
else
    IBS = 2;
end
x(2) = P(2)/K*(1+A*(abs(P(2)/K))^(R-1));

% Loop through all the valuse of P
for i = 3:n
    switch IBS
        
        case 1
            % Reversal Rule
            if (P(i) < P(i-1))
                x_re = x(i-1);  P_re = P(i-1);
                xBC = x_re;     PBC = P_re;
                x(i) = x_re +(P(i)-P_re)/K*(1+A*(abs((P(i)-P_re)/(2*K)))^(R-1));
                if (P(i) > PMAXU)
                    IBS = 5;
                    xMAXU = x_re;
                    PMAXU = P_re;
                    xMAXZ = x_re +(-P_re)/K*(1+A*(abs((-P_re)/(2*K)))^(R-1));
                    xMAXL = -x_re;
                    PMAXL = -P_re;
                else
                    IBS = 6;
                end
            else
                x(i) = P(i)/K*(1+A*(abs(P(i)/K))^(R-1));
            end
            
        case 2 
            % Reversal Rule
            if (P(i) > P(i-1))
                x_re = x(i-1);  P_re = P(i-1);
                xBC = x_re;     PBC = P_re;
                x(i) = x_re +(P(i)-P_re)/K*(1+A*(abs((P(i)-P_re)/(2*K)))^(R-1));
                if (P(i) < PMINL)
                    IBS = 3;
                    xMINL = x_re;
                    PMINL = P_re;
                    xMINZ = x_re +(-P_re)/K*(1+A*(abs((-P_re)/(2*K)))^(R-1));
                    xMINU = -x_re;
                    PMINU = -P_re;
                else
                    IBS = 4;
                end
            else
                x(i) = P(i)/K*(1+A*(abs(P(i)/K))^(R-1));            
            end
            
        case 3
            % Reversal Rule
            if (P(i) < P(i-1))
                x_re = x(i-1);  P_re = P(i-1);
                xZ = x_re +(-P_re)/K*(1+A*(abs((-P_re)/(2*K)))^(R-1));
                if (P(i) > PMAXU)
                    if (xZ > xMAXZ)
                        IBS = 5;
                        xMAXU = x_re;
                        PMAXU = P_re;
                        xMAXZ = xZ;
                        xMAXL = xBC;
                        PMAXL = PBC;
                    else
                        IBS = 11;
                    end
                else
                    IBS = 11;
                end
                xBC = x_re;     PBC = P_re;
            else
                x_re = xMINL;   P_re = PMINL;
            end
            x(i) = x_re +(P(i)-P_re)/K*(1+A*(abs((P(i)-P_re)/(2*K)))^(R-1));
            % Checking intersection with skeletal curve
            xSKEL_prev = P(i-1)/K*(1+A*(abs(P(i-1)/K))^(R-1));
            xSKEL = P(i)/K*(1+A*(abs(P(i)/K))^(R-1));
            if (((xSKEL_prev < x(i-1) && xSKEL >= x(i)) || ...
                (xSKEL_prev > x(i-1) && xSKEL <= x(i))) && ...
                IBS == 3)
                if (P(i) > PMINU)
                    IBS = 1;
                    x(i) = xSKEL;
                else
                    IBS = 3;
                end
            end
            
        case 4
            % Reversal Rule
            if (P(i) < P(i-1))
                x_re = x(i-1);  P_re = P(i-1);
                xZ = x_re +(-P_re)/K*(1+A*(abs((-P_re)/(2*K)))^(R-1));
                if (P(i) > PMAXU)
                    if (xZ > xMAXZ)
                        IBS = 5;
                        xMAXU = x_re;
                        PMAXU = P_re;
                        xMAXZ = xZ;
                        xMAXL = xBC;
                        PMAXL = PBC;                    
                    else
                        IBS = 8;
                    end
                else
                    IBS = 8;
                end
            end
            x(i) = x_re +(P(i)-P_re)/K*(1+A*(abs((P(i)-P_re)/(2*K)))^(R-1));
            % Check intersection with Skeletal curve 
            xSKEL_prev = P(i-1)/K*(1+A*(abs(P(i-1)/K))^(R-1));
            xSKEL = P(i)/K*(1+A*(abs(P(i)/K))^(R-1));
            if (((xSKEL_prev < x(i-1) && xSKEL >= x(i)) || ...
                (xSKEL_prev > x(i-1) && xSKEL <= x(i))) && ...
                IBS == 4)
                if (P(i) > -PBC)
                    IBS = 1;
                    x(i) = xSKEL;
                else
                    IBS = 4;
                end
            end
            
        case 5
            if (P(i) > P(i-1))
                x_re = x(i-1);  P_re = P(i-1);
                xZ = x_re +(-P_re)/K*(1+A*(abs((-P_re)/(2*K)))^(R-1));
                if (P(i) < PMINL)
                    if (xZ < xMINZ)
                        IBS = 3;
                        xMINL = x_re;
                        PMINL = P_re;
                        xMINZ = xZ;
                        xMINU = xBC;
                        PMINU = PBC;
                    else
                        IBS = 12;
                    end
                else
                    IBS = 12;
                end
                xBC = x_re;     PBC = P_re;
            else
                x_re = xMAXU;   P_re = PMAXU;
            end
            x(i) = x_re +(P(i)-P_re)/K*(1+A*(abs((P(i)-P_re)/(2*K)))^(R-1));
            % Intersection Rule for skeletal curve
            xSKEL_prev = P(i-1)/K*(1+A*(abs(P(i-1)/K))^(R-1));
            xSKEL = P(i)/K*(1+A*(abs(P(i)/K))^(R-1));
            if (((xSKEL_prev < x(i-1) && xSKEL >= x(i)) || ...
                (xSKEL_prev > x(i-1) && xSKEL <= x(i))) && ...
                IBS == 5)
                if (P(i) < PMAXL)
                    IBS = 2;
                    x(i) = xSKEL;
                else
                    IBS = 5;
                end
            end
            
        case 6
            % Reversal Rule
            if (P(i) > P(i-1))
                x_re = x(i-1);  P_re = P(i-1);
                xZ = x_re +(-P_re)/K*(1+A*(abs((-P_re)/(2*K)))^(R-1));
                if (P(i) < PMINL)
                    if (xZ < xMINZ)
                        IBS = 3;
                        xMINL = x_re;
                        PMINL = P_re;
                        xMINZ = xZ;
                        xMINU = xBC;
                        PMINU = PBC;
                    else
                        IBS = 7;
                    end
                else
                    IBS = 7;
                end
            end
            x(i) = x_re +(P(i)-P_re)/K*(1+A*(abs((P(i)-P_re)/(2*K)))^(R-1));
            % Intersection Rule for skeletal curve
            xSKEL_prev = P(i-1)/K*(1+A*(abs(P(i-1)/K))^(R-1));
            xSKEL = P(i)/K*(1+A*(abs(P(i)/K))^(R-1));
            if (((xSKEL_prev < x(i-1) && xSKEL >= x(i)) || ...
                (xSKEL_prev > x(i-1) && xSKEL <= x(i))) && ...
                IBS == 6)
                if (P(i) < -PBC)
                    IBS = 2;
                    x(i) = xSKEL;
                else
                    IBS = 6;
                end
            end        
            
        case 7
            % Reversal Rule
            if (P(i) < P(i-1))
                x_re = x(i-1);  P_re = P(i-1);
                IBS = 14;
            end
            x(i) = x_re +(P(i)-P_re)/K*(1+A*(abs((P(i)-P_re)/(2*K)))^(R-1));            
            % Intersection Rule with skeletal curve
            xSKEL_prev = P(i-1)/K*(1+A*(abs(P(i-1)/K))^(R-1));
            xSKEL = P(i)/K*(1+A*(abs(P(i)/K))^(R-1));
            if (((xSKEL_prev < x(i-1) && xSKEL >= x(i)) || ...
                (xSKEL_prev > x(i-1) && xSKEL <= x(i))) && ...
                IBS == 7)
                if (P(i) > PBC)
                    IBS = 1;
                    x(i) = xSKEL;
                else
                    IBS = 7;
                end
            end
            
        case 8
            % Reversal Rule
            if (P(i) > P(i-1))
                x_re = x(i-1);  P_re = P(i-1);
                IBS = 13;
            end
            x(i) = x_re +(P(i)-P_re)/K*(1+A*(abs((P(i)-P_re)/(2*K)))^(R-1));            
            % Intersection Rule with skeletal curve
            xSKEL_prev = P(i-1)/K*(1+A*(abs(P(i-1)/K))^(R-1));
            xSKEL = P(i)/K*(1+A*(abs(P(i)/K))^(R-1));
            if (((xSKEL_prev < x(i-1) && xSKEL >= x(i)) || ...
                (xSKEL_prev > x(i-1) && xSKEL <= x(i))) && ...
                IBS == 8)
                if (P(i) < PBC)
                    IBS = 2;
                    x(i) = xSKEL;
                else
                    IBS = 8;
                end
            end
            
        case 9
            % Reversal Rule
            if (P(i) < P(i-1))
                x_re = x(i-1);  P_re = P(i-1);
                IBS = 14;
            end
            x(i) = x_re +(P(i)-P_re)/K*(1+A*(abs((P(i)-P_re)/(2*K)))^(R-1));            
            % Intersection Rule with min bc
            x3_prev = xMINL +(P(i-1)-PMINL)/K*(1+A*(abs((P(i-1)-PMINL)/(2*K)))^(R-1));
            x3 = xMINL +(P(i)-PMINL)/K*(1+A*(abs((P(i)-PMINL)/(2*K)))^(R-1));
            if (((x3_prev < x(i-1) && x3 >= x(i)) || ...
                (x3_prev > x(i-1) && x3 <= x(i))) && ...
                IBS == 9)
                if (P(i) > PBC)
                    IBS = 3;
                    x(i) = x3;
                else
                    IBS = 9;
                end
            end
            
        case 10
            % Reversal Rule
            if (P(i) > P(i-1))
                x_re = x(i-1);  P_re = P(i-1);
                IBS = 13;
            end
            x(i) = x_re +(P(i)-P_re)/K*(1+A*(abs((P(i)-P_re)/(2*K)))^(R-1));            
            % Intersection Rule with max bc
            x5_prev = xMAXU +(P(i-1)-PMAXU)/K*(1+A*(abs((P(i-1)-PMAXU)/(2*K)))^(R-1));
            x5 = xMAXU +(P(i)-PMAXU)/K*(1+A*(abs((P(i)-PMAXU)/(2*K)))^(R-1));
            if (((x5_prev < x(i-1) && x5 >= x(i)) || ...
                (x5_prev > x(i-1) && x5 <= x(i))) && ...
                IBS == 10)
                if (P(i) < PBC)
                    IBS = 5;
                    x(i) = x5;
                else
                    IBS = 10;
                end
            end
            
        case 11
            % Reversal Rule
            if (P(i) > P(i-1))
                x_re = x(i-1);  P_re = P(i-1);
                IBS = 9;
            end
            x(i) = x_re +(P(i)-P_re)/K*(1+A*(abs((P(i)-P_re)/(2*K)))^(R-1));            
            % Intersection Rule
            xSKEL_prev = P(i-1)/K*(1+A*(abs(P(i-1)/K))^(R-1));            
            xSKEL = P(i)/K*(1+A*(abs(P(i)/K))^(R-1));
            x5_prev = xMAXU +(P(i-1)-PMAXU)/K*(1+A*(abs((P(i-1)-PMAXU)/(2*K)))^(R-1));
            x5 = xMAXU +(P(i)-PMAXU)/K*(1+A*(abs((P(i)-PMAXU)/(2*K)))^(R-1));
            if (((xSKEL_prev < x(i-1) && xSKEL >= x(i)) || ...
                (xSKEL_prev > x(i-1) && xSKEL <= x(i)) || ...
                (x5_prev < x(i-1) && x5 >= x(i)) || ...
                (x5_prev > x(i-1) && x5 <= x(i))) && ...
                IBS == 11)
                if (P(i) < PMINL)
                   if (PMINL < PMAXL)
                       IBS = 2;
                       x(i) = xSKEL;
                   else
                       if (x(i) > x5)
                           IBS = 5;
                           x(i) = x5;
                       else
                           INS = 11;
                       end
                   end
                else
                    IBS = 11;
                end
            end
            
        case 12
            % Reversal Rule
            if (P(i) < P(i-1))
                x_re = x(i-1);  P_re = P(i-1);
                IBS = 10;
            end
            x(i) = x_re +(P(i)-P_re)/K*(1+A*(abs((P(i)-P_re)/(2*K)))^(R-1));
            % Intersection Rule with skeletal curve
            xSKEL_prev = P(i-1)/K*(1+A*(abs(P(i-1)/K))^(R-1)); 
            xSKEL = P(i)/K*(1+A*(abs(P(i)/K))^(R-1));
            x3_prev = xMINL +(P(i-1)-PMINL)/K*(1+A*(abs((P(i-1)-PMINL)/(2*K)))^(R-1));
            x3 = xMINL +(P(i)-PMINL)/K*(1+A*(abs((P(i)-PMINL)/(2*K)))^(R-1));
            if (((xSKEL_prev < x(i-1) && xSKEL >= x(i)) || ...
                (xSKEL_prev > x(i-1) && xSKEL <= x(i)) || ...
                (x3_prev < x(i-1) && x3 >= x(i)) || ...
                (x3_prev > x(i-1) && x3 <= x(i))) && ...
                IBS == 12)
                if (P(i) > PMAXU)
                    if (PMAXU > PMINU)
                        IBS = 1;
                        x(i) = xSKEL;
                    else
                        if (x < x3)
                            IBS = 3;
                            x(i) = x3;
                        else
                            IBS = 12;
                        end
                    end
                else
                    IBS = 12;
                end
            end
            
        case 13
            % Reversal Rule
            if (P(i) < P(i-1))
                x_re = x(i-1);  P_re = P(i-1);
                xZ = x_re +(-P_re)/K*(1+A*(abs((-P_re)/(2*K)))^(R-1));
                if (P(i) > PMAXU)
                    if (xZ > xMAXZ)
                        IBS = 5;
                        xMAXU = x_re;
                        PMAXU = P_re;
                        xMAXZ = xZ;
                        xMAXL = xBC;
                        PMAXL = PBC;  
                    else
                        IBS = 14;
                    end
                else
                    IBS = 14;
                end
            end
            x(i) = x_re +(P(i)-P_re)/K*(1+A*(abs((P(i)-P_re)/(2*K)))^(R-1));
            % Intersection Rule with skeletal curve
            xSKEL_prev = P(i-1)/K*(1+A*(abs(P(i-1)/K))^(R-1)); 
            xSKEL = P(i)/K*(1+A*(abs(P(i)/K))^(R-1));
            x3_prev = xMINL +(P(i-1)-PMINL)/K*(1+A*(abs((P(i-1)-PMINL)/(2*K)))^(R-1));
            x3 = xMINL +(P(i)-PMINL)/K*(1+A*(abs((P(i)-PMINL)/(2*K)))^(R-1));
            if (((xSKEL_prev < x(i-1) && xSKEL >= x(i)) || ...
                (xSKEL_prev > x(i-1) && xSKEL <= x(i)) || ...
                (x3_prev < x(i-1) && x3 >= x(i)) || ...
                (x3_prev > x(i-1) && x3 <= x(i))) && ...
                IBS ==13)
                if (P(i) < PMINU)
                    if (x(i) < x3)
                        IBS = 3;
                        x(i) = x3;
                    else
                        IBS = 13;
                    end
                else
                    if (x(i) < xSKEL)
                        IBS = 1;
                        x(i) = xSKEL;
                    else
                        IBS = 13;
                    end
                end
            end
            
        case 14
            % Reversal Rule
            if (P(i) > P(i-1))
                x_re = x(i-1);  P_re = P(i-1);
                xZ = x_re +(-P_re)/K*(1+A*(abs((-P_re)/(2*K)))^(R-1));
                if (P(i) < PMINL)
                    if (xZ < xMINZ)
                        IBS = 3;
                        xMINL = x_re;
                        PMINL = P_re;
                        xMINZ = xZ;
                        xMINU = xBC;
                        PMINU = PBC;
                    else
                        IBS = 13;
                    end
                else
                    IBS = 13;
                end
            end
            x(i) = x_re +(P(i)-P_re)/K*(1+A*(abs((P(i)-P_re)/(2*K)))^(R-1));
            % Intersection Rule with skeletal curve
            xSKEL_prev = P(i-1)/K*(1+A*(abs(P(i-1)/K))^(R-1));            
            xSKEL = P(i)/K*(1+A*(abs(P(i)/K))^(R-1));
            x5_prev = xMAXU +(P(i-1)-PMAXU)/K*(1+A*(abs((P(i-1)-PMAXU)/(2*K)))^(R-1));
            x5 = xMAXU +(P(i)-PMAXU)/K*(1+A*(abs((P(i)-PMAXU)/(2*K)))^(R-1));
            if (((xSKEL_prev < x(i-1) && xSKEL >= x(i)) || ...
                (xSKEL_prev > x(i-1) && xSKEL <= x(i)) || ...
                (x5_prev < x(i-1) && x5 >= x(i)) || ...
                (x5_prev > x(i-1) && x5 <= x(i))) && ...
                IBS == 14)
                if (P(i) > PMAXL)
                    if (x(i) > x5)
                        IBS = 5;
                        x(i) = x5;
                    else
                        IBS = 14;
                    end
                else
                    if (x(i) > xSKEL)
                        IBS = 2;
                        x(i) = xSKEL;
                    else
                        IBS = 14;
                    end
                end
            end
                
    end
    %IBS_ix(i) = IBS;
end

%strain = x(i)
strain = x;

