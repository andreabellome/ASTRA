function [X4, FVAL4] = findPeriAltitude_IV(rrIN, vvIN, plIN, tIN, gam, plT, hmin, hmax, planetsList)

hh = linspace(hmin, hmax);
DM = zeros(length(hh), 1);

for indh = 1:length(hh)   
    [DM(indh,1)] = Intersection_IV(hh(indh), gam, rrIN, vvIN, plIN, tIN, plT, planetsList);
end

[minDM, row] = min(DM);

if minDM < 1e60
    
    if row == 1 || row == 100
        
        X4    = hh(row);
        FVAL4 = DM(row,1);
        
    else
        
        [X4, FVAL4] = fminbnd(@(H) Intersection_IV(H, gam, rrIN, vvIN, plIN, tIN, plT, planetsList), hh(row - 1), hh(row + 1));

    end
    
else
    X4    = NaN;
    FVAL4 = 1e99;
end