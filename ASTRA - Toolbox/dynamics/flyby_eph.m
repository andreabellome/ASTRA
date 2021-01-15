function [rrOU, vvOU] = flyby_eph(rrIN, vvIN, plIN, tIN, h, gam)

[~, vvPL, ~, ~, radPL, muPL] = ephj2000(plIN, tIN);

vvInfIN = vvIN - vvPL;
n_r     = cross(vvInfIN,vvPL)./norm(cross(vvInfIN,vvPL));

vvInfOU = swingby(vvInfIN, (radPL + h), muPL, gam, n_r);

rrOU = rrIN;
vvOU = vvInfOU' + vvPL;

end