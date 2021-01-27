function [DM, rrIN_2, vvIN_2, tPL] = Intersection_IV(H, gam, rrIN, vvIN, plIN, tIN, plT, planetsList)

% this function computes the anomaly difference between the SC and the planet at the fourth intersection
% moreover it computes the SC state and time at the intersection
%
% INPUT : 
% - H    : flyby altitude (km)
% - gam  : flyby plane inclination (rad)
% - rrIN : SC position vector at the incoming leg of the flyby (km)
% - vvIN : SC velocity vector at the incoming leg of the flyby (km/s)
% - plIN : planet to flyby
% - tIN  : epoch of the planetary encounter (MJD2000)
% - plT  : next planet to flyby
% - planetsList : list of available target planets
% 
% OUTPUT : 
% - DM     : anomaly difference between SC and plT at the intersection (rad)
% - rrIN_2 : SC position vector at the incoming leg at plT (km)
% - vvIN_2 : SC velocity vector at the incoming leg at plT (km)
% - tPL    : epoch of the planetary encounter with plT (MJD2000)

mu = 132724487690; % gravitational parameter of the Sun

% perform the flyby
[rrOU, vvOU] = flyby_ciccopl(rrIN, vvIN, plIN, H, gam, mu);

% check you still arrive at PL
pltsT        = generateTargetBodies(planetsList, rrOU, vvOU);
 
if isempty(find(pltsT == plT,1)) % you do not arrive at plT
    
    DM     = 1e99;
    tPL    = NaN;
    rrIN_2 = NaN.*ones(1,3);
    vvIN_2 = NaN.*ones(1,3);
    
else
    
    kep = car2kep([rrOU, vvOU], mu); % keplerian elements of the outgoing transfer
    
    if kep(2) >= 1 % you are on hyperbolic orbit
        
        DM     = 1e99;
        tPL    = NaN;
        rrIN_2 = NaN.*ones(1,3);
        vvIN_2 = NaN.*ones(1,3);
        
    else % you are on elliptical orbit
        
        [TOFs] = timeofflight(kep, plIN, plT);  % compute the time of flights towards four intersections
        tof_4  = TOFs(4); % extract the time of flight towards fourth intersection
        
        [rrIN_2, vvIN_2] = FGKepler_dt(kep, tof_4, mu);  % compute SC state at the intersection
        
        tPL = tIN + tof_4/86400;  % compute the epoch at the intersection
        
        [~, ~, ~, MPL] = ephj2000(plT, tPL); % compute the planet anomaly at the intersection
        [MSC]          = MeanAn(rrIN_2); % compute the SC anomaly at the intersection
        [DM]           = meanAnDifference(MSC, MPL); % compute the anomaly difference at the intersection
        
    end
    
end

if (tPL - tIN) < 50 % the time of flight is too low
    DM     = 1e99;
    tPL    = NaN;
    rrIN_2 = NaN.*ones(1,3);
    vvIN_2 = NaN.*ones(1,3);
end

end
