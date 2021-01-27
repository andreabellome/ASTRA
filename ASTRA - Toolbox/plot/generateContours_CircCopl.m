function [rpscCONT, EscCONT] = generateContours_CircCopl(pl, vInf)

mu = 132724487690;

alpha = deg2rad(linspace(0,180));

% planet properties
[r, ~, ~] = astroConstantsj2000(pl);
vPL       = sqrt(mu/r);

rpscCONT = zeros(length(alpha), 1);
EscCONT  = zeros(length(alpha), 1);
for indi = 1:length(alpha)
    [rpscCONT(indi,:), EscCONT(indi,:)] = SCorbit(alpha(indi),vInf,vPL,r);
end

% if vInf is too high
[~, row] = min(rpscCONT);

rpscCONT = rpscCONT(1:row);
EscCONT  = EscCONT(1:row);

end