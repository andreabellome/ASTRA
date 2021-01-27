function [delta] = rpip2delta(hpip,rPL,muPL,vInf)

% this function computes the turning angle delta with a given flyby
%
% INPUT : 
% -
% - rPL  : equatorial radius of the planet (km)
% - muPL : gravitational parameter of the planet (km3/s2)
% - vInf : infinity velocity at the planet encounter (km/s)
%
% OUTPUT : 
% - delta : turning angle (rad)

rpip  = hpip + rPL;             % pericentre at the hyperbolic passage [km]
eip   = 1 + (rpip*vInf^2)/muPL; % eccentricity at the hyperbolic passage
delta = 2*asin(1/eip);          % turning angle [rad]

end
