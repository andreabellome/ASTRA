function [X3, FVAL3] = findPeriAltitude_III(rrIN, vvIN, plIN, tIN, gam, plT, hmin, hmax, planetsList)

hh = linspace(hmin, hmax);
DM = zeros(length(hh), 1);

for indh = 1:length(hh)   
    [DM(indh,1)] = Intersection_III(hh(indh), gam, rrIN, vvIN, plIN, tIN, plT, planetsList);
end

[minDM, row] = min(DM);

if minDM < 1e60
    
    if row == 1 || row == 100
        
        X3    = hh(row);
        FVAL3 = DM(row,1);
        
    else
        
        [X3, FVAL3] = fminbnd(@(H) Intersection_III(H, gam, rrIN, vvIN, plIN, tIN, plT, planetsList), hh(row - 1), hh(row + 1));

    end
    
else
    X3    = NaN;
    FVAL3 = 1e99;
end

end