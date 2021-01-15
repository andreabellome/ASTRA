function [DM] = meanAnDifference(M1, M2)

if abs(M1 - M2) > pi
    DM = 2*pi - abs(M1 - M2);
else
    DM = abs(M1 - M2);
end
    
end