function plotArrow2TargetOrbit(path, kept)

mu = 132724487690;
AU = 149597870.7;

% perform the last flyby

hold on;

rrIN   = path(end,1:3);
vvIN   = path(end,4:6);
plIN   = path(end,7);
gamma  = path(end,11);
hp     = path(end,12);
[rrOU, vvOU] = flyby_ciccopl(rrIN, vvIN, plIN, hp, gamma, mu);
kep          = car2kep([rrOU, vvOU], mu);

arrow([kep(1)*(1 - kep(2))/AU, -mu/(2*kep(1))], [kept(1)*(1 - kept(2))/AU, -mu/(2*kept(1))]);

end