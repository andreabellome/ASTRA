function [r, radPL, muPL] = astroConstantsj2000(pl)

% this function computes physiscal properties of a given planet
% 
% INPUT :
% - pl : planet ID (1 - Mercury, 2 - Venus, 3 - Earth, 4 - Mars, 5 - Jupiter, 6 - Saturn, 7 - Uranus, 8 - Neptune)
%
% OUTPUT :
% - r     : planet heliocentric radius (km)
% - radPL : radius of the planet (km)
% - muPL  : gravitational parameter of the planet (km3/s2)
% -------------------------------------------------------------------------

% constants of motion
AU = 149597870.7; % astronomical unit

% initial mean anomaly for the planets
R     = [0.3871*AU 0.7233*AU 1.0000*AU 1.5237*AU 5.2026*AU 9.5549*AU 19.2185*AU 30.1104*AU];

% physical constants for the planets
MUPL  = [22033 3.2486e+05 3.986e+05 42828 1.2669e+08 3.7931e+07 5.7939e+06 6.8347e+06]; % gravitational parameters
RADPL = [2440 6051.8 6378.2 3389.9 69911 58232 25362 24624];                            % radius

r     = R(pl);     % heliocentric circular radius for the given planet
muPL  = MUPL(pl);  % gravitational parameter for the given planet
radPL = RADPL(pl); % radius for the given planet

end
