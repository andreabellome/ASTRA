function [X1, FVAL1] = findPeriAltitude_I(rrIN, vvIN, plIN, tIN, gam, plT, hmin, hmax, planetsList)

hh = linspace(hmin, hmax);
DM = zeros(length(hh), 1);

for indh = 1:length(hh)   
    [DM(indh,1)] = Intersection_I(hh(indh), gam, rrIN, vvIN, plIN, tIN, plT, planetsList);
end

[minDM, row] = min(DM);

if minDM < 1e60
    
    if row == 1 || row == 100
        
        X1    = hh(row);
        FVAL1 = DM(row,1);
        
    else
        
        [X1, FVAL1] = fminbnd(@(H) Intersection_I(H, gam, rrIN, vvIN, plIN, tIN, plT, planetsList), hh(row - 1), hh(row + 1));

    end
    
else
    X1    = NaN;
    FVAL1 = 1e99;
end

end