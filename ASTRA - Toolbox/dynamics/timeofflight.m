function [TOFs] = timeofflight(kep, id_uno, id_due)

% finds the TOF until the next four intersections between two planets

mu    = 132724487690;
r_due = astroConstantsj2000(id_due);

p = kep(1)*(1 - kep(2)^2);
a = kep(1);
e = kep(2);
T = 2*pi*sqrt(a^3/mu);

thSC = kep(6);

TOFs = zeros(1,4);

if id_uno == id_due
    % tempo di volo con lo stesso pianeta
    
    if thSC < pi
        
        % time from periapses passage of the SC position
        MSC = th2M(thSC, e);
        ESC = Mean2E(MSC,e);
        tSC = T/(2*pi)*(ESC - e*sin(ESC));
        
        TOFs(1) = T - 2*tSC;

    else
        
        % time from periapses passage of the SC position
        MSC = th2M(thSC, e);
        ESC = Mean2E(MSC,e);
        tSC = T/(2*pi)*(ESC - e*sin(ESC));
        
        tSC = abs(tSC);
        
        TOFs(1) = 2*tSC;
        
    end
    
    TOFs(2) = T;
    
else

% valuta se stai andando giu o su

if id_uno > id_due
    % stai andando giu
    
    if thSC < pi
        
        % time from periapses passage of the SC position
        MSC = th2M(thSC, e);
        ESC = Mean2E(MSC,e);
        tSC = T/(2*pi)*(ESC - e*sin(ESC));

        % true anomaly at the planet encounter
        argth   = (p - r_due)/(e*r_due);
        [argth] = checkArgTheta(argth);
        thPL    = acos(argth);
        
        % time from periapses passage of the encounter position
        MPL = th2M(thPL, e);
        EPL = Mean2E(MPL,e);
        tPL = T/(2*pi)*(EPL - e*sin(EPL));

        TOFs(1) = T - tSC - tPL;
        TOFs(2) = TOFs(1) + 2*tPL;
        
    else
        
        % time from periapses passage of the SC position
        MSC = th2M(thSC, e);
        ESC = Mean2E(MSC,e);
        tSC = T/(2*pi)*(ESC - e*sin(ESC));
        
        % true anomaly at the planet encounter
        argth   = (p - r_due)/(e*r_due);
        [argth] = checkArgTheta(argth);
        thPL    = acos(argth);
        
        % time from periapses passage of the encounter position
        MPL = th2M(thPL, e);
        EPL = Mean2E(MPL,e);
        tPL = T/(2*pi)*(EPL - e*sin(EPL));
        
        TOFs(1) = abs(tSC) - tPL;
        TOFs(2) = TOFs(1) + 2*tPL;
        
    end
    
else
    % stai andando su
    
    if thSC >= 0 && thSC < pi
        
        % time from periapses passage of the SC position
        MSC = th2M(thSC, e);
        ESC = Mean2E(MSC,e);
        tSC = T/(2*pi)*(ESC - e*sin(ESC));
        
        % true anomaly at the planet encounter
        argth   = (p - r_due)/(e*r_due);
        [argth] = checkArgTheta(argth);
        thPL    = acos(argth);
        
        % time from periapses passage of the encounter position
        MPL = th2M(thPL, e);
        EPL = Mean2E(MPL,e);
        tPL = T/(2*pi)*(EPL - e*sin(EPL));
        
        TOFs(1) = tPL - tSC;
        TOFs(2) = T - tSC - tPL;
        
    else
        
        % time from periapses passage of the SC position
        MSC = th2M(thSC, e);
        ESC = Mean2E(MSC,e);
        tSC = T/(2*pi)*(ESC - e*sin(ESC));
        
        % true anomaly at the planet encounter
        argth   = (p - r_due)/(e*r_due);
        [argth] = checkArgTheta(argth);
        thPL    = acos(argth);
        
        % time from periapses passage of the encounter position
        MPL = th2M(thPL, e);
        EPL = Mean2E(MPL,e);
        tPL = T/(2*pi)*(EPL - e*sin(EPL));
        
        tSC = abs(tSC);
        TOFs(1) = tSC + tPL;
        TOFs(2) = T - tPL + tSC;

    end
    
end

end

TOFs(3) = TOFs(1) + T;
TOFs(4) = TOFs(2) + T;

function M = th2M(theta, e)

if e < 1
E = 2*atan(sqrt((1-e)/(1+e))*tan(theta/2));
M =  E - e*sin(E);
elseif e > 1
E = 2*atanh(sqrt((e-1)/(e+1))*tan(theta/2));
M = e*sinh(E) - E;
end

end

function E=Mean2E(M,e)
%
tol   = 1e-13;
err   = 1;
E     = M+e*cos(M);   %initial guess
maxit = 1000;
it    = 0;
while (err>tol) && (it<maxit)
    Enew = E-(E-e*sin(E)-M)/(1-e*cos(E));
    err  = abs(E-Enew);
    E    = Enew;
    it   = it + 1;
end

end

function [argth] = checkArgTheta(argth)

if abs(argth - 1) < 1e-2
    argth = 1;
elseif abs(argth + 1) < 1e-2
    argth = -1;
elseif abs(argth) < 1e-2
    argth = 0;
else
    
end

end

end