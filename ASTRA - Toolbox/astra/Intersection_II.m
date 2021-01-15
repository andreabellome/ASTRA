function [DM, rrIN_2, vvIN_2, tPL] = Intersection_II(H, gam, rrIN, vvIN, plIN, tIN, plT, planetsList)

mu = 132724487690;

% perform the flyby
[rrOU, vvOU] = flyby_ciccopl(rrIN, vvIN, plIN, H, gam, mu);

% check you still arrive at PL
% planetsList  = [2 3 4 5];
pltsT        = generateTargetBodies(planetsList, rrOU, vvOU);

if isempty(find(pltsT == plT,1))
    
    DM     = 1e99;
    tPL    = NaN;
    rrIN_2 = NaN.*ones(1,3);
    vvIN_2 = NaN.*ones(1,3);
    
else
    
    kep = car2kep([rrOU, vvOU], mu);
    
    if kep(2) >= 1
        
        DM     = 1e99;
        tPL    = NaN;
        rrIN_2 = NaN.*ones(1,3);
        vvIN_2 = NaN.*ones(1,3);
        
    else
        
        [TOFs] = timeofflight(kep, plIN, plT);
        tof_2  = TOFs(2);
        
        [rrIN_2, vvIN_2] = FGKepler_dt(kep, tof_2, mu);
        
        tPL = tIN + tof_2/86400;
        
        [~, ~, ~, MPL] = ephj2000(plT, tPL);
        [MSC]          = MeanAn(rrIN_2);
        [DM]           = meanAnDifference(MSC, MPL);
        
    end
    
end

if (tPL - tIN) < 50
    DM     = 1e99;
    tPL    = NaN;
    rrIN_2 = NaN.*ones(1,3);
    vvIN_2 = NaN.*ones(1,3);
end

end