function [dv] = find_dvToTargetFromMMAT(MAT_2, kept)

mu = 132724487690;

% find keplerian elements
% [kep] = find_finalKEPsFromMMAT(MAT_2);
[kep] = finalKEPsFromMMAT(MAT_2, mu);

% find dv from target orbit
dv = zeros(size(kep,1), 1);
for indk = 1:size(kep,1)
    dv(indk,1) = DVcost(kep(indk,:), kept);
end

end