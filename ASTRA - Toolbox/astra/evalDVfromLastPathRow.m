function [dv, Tt, POS, DV, ALIGN] = evalDVfromLastPathRow(path, kept)

mu = 132724487690;

rrIN  = path(end,1:3);
vvIN  = path(end,4:6);
plIN  = path(end,7);
gamma = path(end,11);
hp    = path(end,12);

[rrOU, vvOU] = flyby_ciccopl(rrIN, vvIN, plIN, hp, gamma, mu);
kep          = car2kep([rrOU, vvOU], mu);

% [dv] = DVcost_v2(kep, kept);
[dv, Tt, POS, DV, ALIGN] = DVcost_v3(kep, kept);

end