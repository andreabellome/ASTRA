function [rrOU, vvOU] = flyby_ciccopl(rrIN, vvIN, plIN, h, gam, mu)

[rPL, radPL, muPL] = astroConstantsj2000(plIN);

vvPL    = sqrt(mu/rPL).*[-sin(MeanAn(rrIN)) cos(MeanAn(rrIN)) 0];
vvInfIN = vvIN - vvPL;
n_r     = cross(vvInfIN,vvPL)./norm(cross(vvInfIN,vvPL));

vvInfOU = swingby(vvInfIN, (radPL + h), muPL, gam, n_r);

rrOU = rrIN;
vvOU = vvInfOU' + vvPL;

end