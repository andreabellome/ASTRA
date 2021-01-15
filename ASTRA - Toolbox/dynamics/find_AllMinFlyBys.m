function [ascFB, EscFB, rpscFB, TscFB, delta] = find_AllMinFlyBys(pl, vInf)

mu = 132724487690;

% planet properties
[r, rPL, muPL] = astroConstantsj2000(pl);

% minimum altitude flyby of the planet
[hmin]  = maxmin_flybyAltitude(pl);
[delta] = (rpip2delta(hmin, rPL, muPL, vInf));

alpha(1) = 0;
indi     = 2;
while alpha(end) < pi
    alpha(indi) = alpha(indi - 1) + delta;
    indi = indi + 1;
end
alpha(end) = [];

for indi = 1:length(alpha)
%     [ascFB(indi), EscFB(indi), rpscFB(indi), TscFB(indi)] = ...
%         generateOrbit(pl, vInf, alpha(indi));   
    [rpscFB(indi), EscFB(indi)] = SCorbit(alpha(indi), vInf, sqrt(mu/r), r);
    ascFB(indi) = -mu/(2*EscFB(indi));
    TscFB(indi) = 2*pi*sqrt(ascFB(indi)^3/mu);
end

% if vInf is too high
[~, row] = min(rpscFB);

rpscFB = rpscFB(1:row);
EscFB  = EscFB(1:row);
ascFB  = ascFB(1:row);
TscFB  = TscFB(1:row);

end