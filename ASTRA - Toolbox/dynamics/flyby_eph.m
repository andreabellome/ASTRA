function [rrOU, vvOU] = flyby_eph(rrIN, vvIN, plIN, tIN, h, gam)

% this function computes the flyby with a planet given its position in ephemerides model 
% 
% INPUT : 
% - rrIN : position vector of the incoming leg at the planet (km)
% - vvIN : velocity vector of the incoming leg at the planet (km/s)
% - plIN : planet to flyby
% - h    : flyby altitude (km)
% - gam  : flyby hyperbola plane inclination (rad)
% - mu   : gravitational parameter of the Sun (km3/s2)
% 
% OUTPUT : 
% - rrOU : position vector of the outgoing leg at the planet (km)
% - vvOU : velocity vector of the incoming leg at the planet (km/s)

% TBR (Andrea) : check gamma angle
% TBR (Andrea) : use this function or the circular coplanar one?

[~, vvPL, ~, ~, radPL, muPL] = ephj2000(plIN, tIN); % find planet properties from ephemerides

vvInfIN = vvIN - vvPL; % infinity velocity at planet encounter
n_r     = cross(vvInfIN,vvPL)./norm(cross(vvInfIN,vvPL)); % reference vector for the flyby plane

vvInfOU = swingby(vvInfIN, (radPL + h), muPL, gam, n_r); % compute the outgoing infinity velocity

% outgoing parameters
rrOU = rrIN; % outgoing position vector
vvOU = vvInfOU' + vvPL; % outgoing velocity vector

end
