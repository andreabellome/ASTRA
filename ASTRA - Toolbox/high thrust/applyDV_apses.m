function [x, y] = applyDV_apses(kep, dvMax)

mu = 132724487690;

% given the central keplerian elements

sma = kep(1);
ecc = kep(2);

ra  = sma*(1 + ecc);
rp  = sma*(1 - ecc);

va  = sqrt((2*mu)/(ra + rp)*(rp/ra));
vp  = sqrt((2*mu)/(ra + rp)*(ra/rp));

vpSOppm = vp + dvMax;
rrppm    = [rp, 0, 0];
vvppm    = [0, vpSOppm, 0];
kepSOppm = car2kep([rrppm, vvppm], mu);
raSOppm  = kepSOppm(1)*(1 + kepSOppm(2));
ESOppm   = -mu/(2*kepSOppm(1));

vpSOppm = vp - dvMax;
rrppm     = [rp, 0, 0];
vvppm     = [0, vpSOppm, 0];
kepSOppm  = car2kep([rrppm, vvppm], mu);
raSOppm_2 = kepSOppm(1)*(1 + kepSOppm(2));
ESOppm_2  = -mu/(2*kepSOppm(1));

vaSOapm = va + dvMax;

rpSOapm_3 = (2*mu/((vaSOapm^2)*ra) - 1)^(-1)*ra;

rrapm     = [-ra, 0, 0];
vvapm     = [0, -vaSOapm, 0];
kepSOapm  = car2kep([rrapm, vvapm], mu);
raSOppm_3 = kepSOapm(1)*(1 + kepSOapm(2));
ESOapm_3  = -mu/(2*kepSOapm(1));

vaSOapm = va - dvMax;

rpSOapm_4 = (2*mu/((vaSOapm^2)*ra) - 1)^(-1)*ra;

rrapm     = [-ra, 0, 0];
vvapm     = [0, -vaSOapm, 0];
kepSOapm  = car2kep([rrapm, vvapm], mu);
raSOppm_4 = kepSOapm(1)*(1 + kepSOapm(2));
ESOapm_4  = -mu/(2*kepSOapm(1));

% save rp and Energy
x = [rpSOapm_4, rp, rpSOapm_3, rp, rp, rpSOapm_4];
y = [ESOapm_4, ESOppm_2, ESOapm_3, ESOppm, ESOppm, ESOapm_4];

end