function [X2, FVAL2] = findPeriAltitude_II(rrIN, vvIN, plIN, tIN, gam, plT, hmin, hmax, planetsList)

hh = linspace(hmin, hmax);
DM = zeros(length(hh), 1);

for indh = 1:length(hh)   
    [DM(indh,1)] = Intersection_II(hh(indh), gam, rrIN, vvIN, plIN, tIN, plT, planetsList);
end

[minDM, row] = min(DM);

if minDM < 1e60
    
    if row == 1 || row == 100
        
        X2    = hh(row);
        FVAL2 = DM(row,1);
        
    else
        
        [X2, FVAL2] = fminbnd(@(H) Intersection_II(H, gam, rrIN, vvIN, plIN, tIN, plT, planetsList), hh(row - 1), hh(row + 1));

    end
    
else
    X2    = NaN;
    FVAL2 = 1e99;
end

end