function [ascFB, EscFB, rpscFB, TscFB] = find_AllMinFlyBys(pl, vInf)

% this function computes the flyby orbits on a Tisserand contour specified by planet and 
% infinity velocity (orbits are spaced by minimum altitude constraint)
%
% INPUT :
% - pl : planet ID (1 - Mercury, 2 - Venus, 3 - Earth, 4 - Mars, 5 - Jupiter, 6 - Saturn, 7 - Uranus, 8 - Neptune)
% - vInf : infinity velocity magnitude at the planet (km/s)
% 
% OUTPUT : 
% - ascFB  : vector with semi-major axis of flyby orbits (km)
% - EscFB  : vector with energy of flyby orbits (km2/s2)
% - rpscFB : vector with perihelion of flyby orbits (km)
% - TscFB  : vector with period of flyby orbits (secs)

mu = 132724487690; % gravitational parameter of the Sun

% planet properties
[r, rPL, muPL] = astroConstantsj2000(pl);

% minimum altitude flyby of the planet
[hmin]  = maxmin_flybyAltitude(pl);
[delta] = (rpip2delta(hmin, rPL, muPL, vInf)); % turning angle (minimum altitude constraint)

% find alpha angles (spaced by minimum altitude constraint)
alpha(1) = 0;
indi     = 2;
while alpha(end) < pi
    alpha(indi) = alpha(indi - 1) + delta;
    indi = indi + 1;
end
alpha(end) = [];

% for each alpha, find the corresponding spacecraft orbit
for indi = 1:length(alpha)
    [rpscFB(indi), EscFB(indi)] = SCorbit(alpha(indi), vInf, sqrt(mu/r), r); % find spacecraft orbit
    ascFB(indi) = -mu/(2*EscFB(indi)); % semi-major axis
    TscFB(indi) = 2*pi*sqrt(ascFB(indi)^3/mu); % period
end

% if vInf is too high
[~, row] = min(rpscFB);

rpscFB = rpscFB(1:row);
EscFB  = EscFB(1:row);
ascFB  = ascFB(1:row);
TscFB  = TscFB(1:row);

end
