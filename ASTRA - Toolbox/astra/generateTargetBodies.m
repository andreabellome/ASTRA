function [planetsTarget] = generateTargetBodies(planetsList, rr, vv)

mu = 132724487690;

if norm(rr) > 1e98 || ~isreal(rr) || ~isreal(vv)
    planetsTarget = [];
else

    sma = -mu/(2*(norm(vv)^2/2 - mu/norm(rr)));
    ecc = norm(-rr./norm(rr) + 1/mu*cross(vv,cross(rr,vv)));
    rp  = sma*(1 - ecc);
    ra  = sma*(1 + ecc);
    
    planetsTarget = [];
    
    for indi = 1:length(planetsList)
        
        PL  = planetsList(indi);
        rPL = astroConstantsj2000(PL);
        
        if (rp <= rPL && ra >= rPL) || abs(rp - rPL)<1e-6 || abs(ra - rPL)<1e-6
            planetsTarget = [planetsTarget planetsList(indi)];
        end
        
    end
end

end