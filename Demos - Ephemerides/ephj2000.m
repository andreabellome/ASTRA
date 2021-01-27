function [rr, vv, r, M, radPL, muPL] = ephj2000(pl, t)

% circular coplanar ephemerides for the planets in MJD2000

% INPUT : 
% pl - planet ID (Mercury, Venus, Earth, Mars, Jupiter, Saturn, Uranus, Neptune)
% t  - epoch MJD2000 (the reference one is January 1, 2000)

% t is expressed in days from t0 = 0

% OUTPUT :
% rr    : planet position at time t
% vv    : planet velocity at time t
% r     : norm(rr)
% M     : planet angular position at time t (0<= M <= 2*pi)
% radPL : radius of the planet
% muPL  : gravitational parameter of the planet

mu = 132724487690;
AU = 149597870.7;

R     = [0.3871*AU 0.7233*AU 1.0000*AU 1.5237*AU 5.2026*AU 9.5549*AU 19.2185*AU 30.1104*AU];
M0    = [4.40261049144822 3.17615017277928 1.75346248630862 6.203476120241 0.599538051352572 0.8740085295212 5.48129378235079 5.31189212515222];
MUPL  = [22033 3.2486e+05 3.986e+05 42828 1.2669e+08 3.7931e+07 5.7939e+06 6.8347e+06];
RADPL = [2440 6051.8 6378.2 3389.9 69911 58232 25362 24624];

r     = R(pl);
M0    = M0(pl);
muPL  = MUPL(pl);
radPL = RADPL(pl);

% angular velocity
om = sqrt(mu/r)/r;                   % rad/s

% planet position and velocity at time t
M  = M0 + om*(t*86400);              % rad
M  = wrapTo2Pi(M);                   % rad
rr = r.*[cos(M) sin(M) 0];           % km
vv = sqrt(mu/r).*[-sin(M) cos(M) 0]; % km/s

end