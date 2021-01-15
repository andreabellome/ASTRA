function [delta] = rpip2delta(hpip,rPL,muPL,vInf)

rpip  = hpip + rPL;             % pericentre at the hyperbolic passage [km]
eip   = 1 + (rpip*vInf^2)/muPL; % eccentricity at the hyperbolic passage
delta = 2*asin(1/eip);          % turning angle at the hyperbolic passage [rad]

end