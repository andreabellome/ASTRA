function [dv, Tt, POS, DV, ALIGN] = evalDVfromLastPathRow(path, kept)

% this function computes the Hohmann DV to arrive to target orbit from last planetary encounter
%
% INPUT :
% - path : matrix containing all the info for a given planetary sequence
%          - row 1 and cols 1:6     -> initial state at the departure
%          - row 2:end and cols 1:6 -> incoming state at the planet encounter
%          - row 1:end and cols 7   -> planetary sequence
%          - row 1:end and cols 8   -> epochs of planetary encounters
%          - row 2:end and cols 11  -> flyby inclinations
%          - row 2:end and cols 12  -> flyby altitudes
% - kept : target orbit keplerian elements
% 
% OUTPUT :
% - dv    : total DV (km/s)
% - Tt    : Hohmann transfer time (secs)
% - POS   : position of Hohmann manoeuvres (periapsis or apoapsis)
% - DV    : magnitude of Hohmann manoeuvres (km/s)
% - ALIGN : alignment of Hohmann manoeuvres (parallel or antiparallel the velocity)

mu = 132724487690; % gravitational parameter of the Sun

rrIN  = path(end,1:3); % SC position vector for incoming leg (km)
vvIN  = path(end,4:6); % SC velocity vector for incoming leg (km)
plIN  = path(end,7);   % planet to flyby
gamma = path(end,11);  % flyby inclination
hp    = path(end,12);  % flyby altitude 

[rrOU, vvOU] = flyby_ciccopl(rrIN, vvIN, plIN, hp, gamma, mu); % compute the flyby with the last planet of the sequence
kep          = car2kep([rrOU, vvOU], mu); % keplerian elements of the outgoing leg

% [dv] = DVcost_v2(kep, kept);
[dv, Tt, POS, DV, ALIGN] = DVcost_v3(kep, kept); % compute the Hohmann DV from kep to target orbit kept 

end
