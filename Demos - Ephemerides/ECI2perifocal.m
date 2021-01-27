function [posPERI, velPERI, vr, vth, k] = ECI2perifocal(posECI, velECI, mu)

kep  = car2kep([posECI, velECI], mu);
a    = kep(1);
e    = kep(2);
th   = kep(6);

r   = a*(1 - e^2)/(1 + e*cos(th));
vr  = sqrt(mu/(a*(1 - e^2)))*e*sin(th);       % radial velocity
vth = sqrt(mu/(a*(1 - e^2)))*(1 + e*cos(th)); % transverse velocity

posPERI = r.*[cos(th) sin(th) 0];

vvr     = vr.*[cos(th) sin(th) 0];
vvth    = vth.*[-sin(th) cos(th) 0];

velPERI = vvr + vvth;

% compute flight-path angle
k = atan(vr/vth); % defined in [-90 90]

end