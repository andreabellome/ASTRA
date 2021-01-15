function [rrOU, vvOU] = flyby_ciccopl(rrIN, vvIN, plIN, h, gam, mu)

% this function computes the circular-coplanar flyby
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
% - vvOU : velocity vector of the incoming leg at the planet (km)

% TBR (Andrea) : check gamma angle

[rPL, radPL, muPL] = astroConstantsj2000(plIN); % physical constants of the given planet

vvPL    = sqrt(mu/rPL).*[-sin(MeanAn(rrIN)) cos(MeanAn(rrIN)) 0]; % planet heliocentric velocity vector
vvInfIN = vvIN - vvPL; % infinity velocity at the planet encounter
n_r     = cross(vvInfIN,vvPL)./norm(cross(vvInfIN,vvPL)); % reference vector for the flyby plane

vvInfOU = swingby(vvInfIN, (radPL + h), muPL, gam, n_r); % compute the outgoing infinity velocity

% outgoing parameters
rrOU = rrIN; % outgoing position vector
vvOU = vvInfOU' + vvPL; % outgoing velocity vector

end
