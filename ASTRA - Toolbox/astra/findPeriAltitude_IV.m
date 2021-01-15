function [X4, FVAL4] = findPeriAltitude_IV(rrIN, vvIN, plIN, tIN, gam, plT, hmin, hmax, planetsList)

% this function computes the flyby altitude at a given planet (plIN) to encounter the next planet (plT)
% at the fourth intersection
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
% - X4    : flyby altitude (km)
% - FVAL4 : anomaly difference between SC and plT at the intersection (rad)

% check for DM=0 -> used to initialise the fminbnd search
hh = linspace(hmin, hmax);
DM = zeros(length(hh), 1);
for indh = 1:length(hh)   
    [DM(indh,1)] = Intersection_IV(hh(indh), gam, rrIN, vvIN, plIN, tIN, plT, planetsList);
end
[minDM, row] = min(DM); % select the minimum one

if minDM < 1e60
    
    if row == 1 || row == 100 % the solution is either the min. or max. altitude
        
        X4    = hh(row);
        FVAL4 = DM(row,1);
        
    else
        
         % solve the phasing problem with fminbnd
        [X4, FVAL4] = fminbnd(@(H) Intersection_IV(H, gam, rrIN, vvIN, plIN, tIN, plT, planetsList), hh(row - 1), hh(row + 1));

    end
    
else % no solution exists
    X4    = NaN;
    FVAL4 = 1e99;
end
