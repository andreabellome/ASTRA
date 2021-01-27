function [TOFs] = timeofflight(kep, id_uno, id_due)

% finds the TOF until the next four intersections between two planets (circular coplanar model)
% 
% INPUT : 
% - kep    : keplerian elements of the transfer orbit [1x6] vector
% - id_uno : first planet ID (1 - Mercury, 2 - Venus, 3 - Earth, 4 - Mars, 5 - Jupiter, 6 - Saturn, 7 - Uranus, 8 - Neptune)
% - id_due : second planet ID (1 - Mercury, 2 - Venus, 3 - Earth, 4 - Mars, 5 - Jupiter, 6 - Saturn, 7 - Uranus, 8 - Neptune)
% 
% OUTPUT : 
% - TOFs : time of flight vector towards four intersections (secs)

mu    = 132724487690; % gravitational parameter of the Sun
r_due = astroConstantsj2000(id_due); % second planet heliocentric radius

p = kep(1)*(1 - kep(2)^2); % semi-latus rectum of the transfer orbit
a = kep(1); % semi-major axis of the transfer orbit
e = kep(2); % eccentricity of the transfer orbit
T = 2*pi*sqrt(a^3/mu); % period of the transfer orbit

thSC = kep(6); % intial spacecraft orbit

TOFs = zeros(1,4); % pre-assign the time of flight vector

if id_uno == id_due
    % time of flight for a same planet-to-planet transfer
    
    if thSC < pi % the SC is either in the first or second quadrant
        
        % time from periapses passage of the SC position
        MSC = th2M(thSC, e); % SC mean anomaly
        ESC = Mean2E(MSC,e); % SC eccentric anomaly
        tSC = T/(2*pi)*(ESC - e*sin(ESC)); % time from periapsis passage
        
        TOFs(1) = T - 2*tSC; % time of flight towards the first intersection

    else % the SC is either in the third or fourth quadrant
        
        % time from periapses passage of the SC position
        MSC = th2M(thSC, e); % SC mean anomaly
        ESC = Mean2E(MSC,e); % SC eccentric anomaly
        tSC = T/(2*pi)*(ESC - e*sin(ESC)); % time from periapsis passage
        tSC = abs(tSC); % non-negative time from periapsis passage
        
        TOFs(1) = 2*tSC; % time of flight towards the first intersection
        
    end
    
    TOFs(2) = T; % time of flight towards the second intersection
    
else

% check if you are going up or down

if id_uno > id_due
    % you are going down (e.g. from Earth to Venus)
    
    if thSC < pi % the SC is either in the first or second quadrant
        
        % time from periapses passage of the SC position
        MSC = th2M(thSC, e); % SC mean anomaly
        ESC = Mean2E(MSC,e); % SC eccentric anomaly
        tSC = T/(2*pi)*(ESC - e*sin(ESC)); % time from periapsis passage

        % true anomaly at the planet encounter
        argth   = (p - r_due)/(e*r_due); % cos(th) = argth
        [argth] = checkArgTheta(argth); % wrap argth within tolerances (0 and pi)
        thPL    = acos(argth); % true anomaly at the planet encounter
        
        % time from periapses passage of the encounter position
        MPL = th2M(thPL, e); % SC mean anomaly at encounter position
        EPL = Mean2E(MPL,e); % SC eccentric anomaly at encounter position
        tPL = T/(2*pi)*(EPL - e*sin(EPL)); % time from periapsis passage of encounter position

        TOFs(1) = T - tSC - tPL; % time of flight towards the first intersection
        TOFs(2) = TOFs(1) + 2*tPL; % time of flight towards the second intersection
        
    else % the SC is either in the third or fourth quadrant
        
        % time from periapses passage of the SC position
        MSC = th2M(thSC, e); % SC mean anomaly
        ESC = Mean2E(MSC,e); % SC eccentric anomaly
        tSC = T/(2*pi)*(ESC - e*sin(ESC)); % time from periapsis passage
        
        % true anomaly at the planet encounter
        argth   = (p - r_due)/(e*r_due); % cos(th) = argth
        [argth] = checkArgTheta(argth); % wrap argth within tolerances (0 and pi)
        thPL    = acos(argth);  % true anomaly at the planet encounter
        
        % time from periapses passage of the encounter position
        MPL = th2M(thPL, e); % SC mean anomaly at encounter position
        EPL = Mean2E(MPL,e); % SC eccentric anomaly at encounter position
        tPL = T/(2*pi)*(EPL - e*sin(EPL)); % time from periapsis passage of encounter position
        
        TOFs(1) = abs(tSC) - tPL; % time of flight towards the first intersection
        TOFs(2) = TOFs(1) + 2*tPL; % time of flight towards the second intersection
        
    end
    
else
    % you are going uo (e.g. from Venus to Earth)
    
    if thSC >= 0 && thSC < pi % the SC is either in the first or second quadrant
        
        % time from periapses passage of the SC position
        MSC = th2M(thSC, e); % SC mean anomaly
        ESC = Mean2E(MSC,e); % SC eccentric anomaly
        tSC = T/(2*pi)*(ESC - e*sin(ESC)); % time from periapsis passage
        
        % true anomaly at the planet encounter
        argth   = (p - r_due)/(e*r_due); % cos(th) = argth
        [argth] = checkArgTheta(argth); % wrap argth within tolerances (0 and pi)
        thPL    = acos(argth); % true anomaly at the planet encounter
        
        % time from periapses passage of the encounter position
        MPL = th2M(thPL, e);  % SC mean anomaly at encounter position
        EPL = Mean2E(MPL,e);  % SC eccentric anomaly at encounter position
        tPL = T/(2*pi)*(EPL - e*sin(EPL)); % time from periapsis passage of the encounter position
        
        TOFs(1) = tPL - tSC; % time of flight towards the first intersection 
        TOFs(2) = T - tSC - tPL; % time of flight towards the second intersection 
        
    else % the SC is either in the third or fourth quadrant
        
        % time from periapses passage of the SC position
        MSC = th2M(thSC, e); % SC mean anomaly
        ESC = Mean2E(MSC,e); % SC eccentric anomaly
        tSC = T/(2*pi)*(ESC - e*sin(ESC)); % time from periapsis passage
        
        % true anomaly at the planet encounter
        argth   = (p - r_due)/(e*r_due); % cos(th) = argth
        [argth] = checkArgTheta(argth); % wrap argth within tolerances (0 and pi)
        thPL    = acos(argth); % true anomaly at the planet encounter
        
        % time from periapses passage of the encounter position
        MPL = th2M(thPL, e); % SC mean anomaly at encounter position
        EPL = Mean2E(MPL,e); % SC eccentric anomaly at encounter position
        tPL = T/(2*pi)*(EPL - e*sin(EPL)); % time from periapsis passage of encounter position
        
        tSC = abs(tSC); % non-negative time from periapsis passage
        TOFs(1) = tSC + tPL; % time of flight towards the first intersection 
        TOFs(2) = T - tPL + tSC; % time of flight towards the second intersection 

    end
    
end

end

TOFs(3) = TOFs(1) + T; % time of flight towards the third intersection 
TOFs(4) = TOFs(2) + T; % time of flight towards the fourth intersection 

% -----------------------------------------------------
% additional functions are included below

function M = th2M(theta, e)

% this function computes the mean anomaly from true anomaly and eccentricity
% 
% INPUT : 
% - theta : true anomaly (rad)
% - e     : eccentricity
% 
% OUTPUT:
% - M : mean anomaly (rad)

if e < 1
E = 2*atan(sqrt((1-e)/(1+e))*tan(theta/2));
M =  E - e*sin(E);
elseif e > 1
E = 2*atanh(sqrt((e-1)/(e+1))*tan(theta/2));
M = e*sinh(E) - E;
end

end

function E=Mean2E(M,e)

% this function computes the eccentric anomaly from mean anomaly and eccentricity
% 
% INPUT : 
% - M : mean anomaly (rad)
% - e : eccentricity
% 
% OUTPUT:
% - E : eccentric anomaly (rad)

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

% this function computes the cos(th)=argth given some tolerances

if abs(argth - 1) < 1e-2 % tolerance on pi
    argth = 1;
elseif abs(argth + 1) < 1e-2 % tolerance on pi
    argth = -1;
elseif abs(argth) < 1e-2 % tolerance on 0 and 2pi
    argth = 0;
else
    
end

end

end
