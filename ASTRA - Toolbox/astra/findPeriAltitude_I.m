function [X1, FVAL1] = findPeriAltitude_I(rrIN, vvIN, plIN, tIN, gam, plT, hmin, hmax, planetsList)

% this function computes the flyby altitude at a given planet (plIN) to encounter the next planet (plT)
% at the first intersection
%
% INPUT :
% - rrIN : SC position vector at the incoming leg of the flyby (km)
% - vvIN : SC velocity vector at the incoming leg of the flyby (km/s)
% - plIN : planet to flyby
% - tIN  : epoch of the planetary encounter (MJD2000)
% - gam  : flyby plane inclination (rad)
% - plT  : next planet to flyby
% - hmin : minimum flyby altitude (km)
% - hmax : maximum flyby altitude (km)
% - planetsList : list of available target planets
%
% OUTPUT :
% - X1    : flyby altitude (km)
% - FVAL1 : anomaly difference between SC and plT at the intersection (rad)

% check for DM=0 -> used to initialise the fminbnd search
hh = linspace(hmin, hmax);
DM = zeros(length(hh), 1);
for indh = 1:length(hh)   
    [DM(indh,1)] = Intersection_I(hh(indh), gam, rrIN, vvIN, plIN, tIN, plT, planetsList);
end
[minDM, row] = min(DM); % select the minimum one

if minDM < 1e60
    
    if row == 1 || row == 100 % the solution is either the min. or max. altitude
        
        X1    = hh(row);
        FVAL1 = DM(row,1);
        
    else
        
        % solve the phasing problem with fminbnd
        [X1, FVAL1] = fminbnd(@(H) Intersection_I(H, gam, rrIN, vvIN, plIN, tIN, plT, planetsList), hh(row - 1), hh(row + 1));

    end
    
else % no solution exists
    X1    = NaN;
    FVAL1 = 1e99;
end

end
