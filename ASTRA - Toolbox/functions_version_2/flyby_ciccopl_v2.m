function [rrOU, vvOU] = flyby_ciccopl_v2(rrIN, vvIN, plIN, h, gam, mu)

% h is a vector

[rPL, radPL, muPL] = astroConstantsj2000(plIN);

vvPL    = sqrt(mu/rPL).*[-sin(MeanAn(rrIN)) cos(MeanAn(rrIN)) 0];
vvInfIN = vvIN - vvPL;

if norm(cross(vvInfIN, vvPL)) < 1e-5
    n_r = [0 0 1];
else
    n_r = cross(vvInfIN, vvPL)./norm(cross(vvInfIN,vvPL));
end
% n_r     = cross(vvInfIN, vvPL)./norm(cross(vvInfIN,vvPL));
vvInfOU = swingby_v2(vvInfIN, (radPL + h), muPL, gam, n_r);

rrOU      = rrIN.*ones(size(vvInfOU,1),1);
vvOU(:,1) = vvInfOU(:,1) + vvPL(1);
vvOU(:,2) = vvInfOU(:,2) + vvPL(2);
vvOU(:,3) = vvInfOU(:,3) + vvPL(3);

end