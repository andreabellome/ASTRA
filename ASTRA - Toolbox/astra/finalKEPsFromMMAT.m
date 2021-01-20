function [kep] = finalKEPsFromMMAT(MAT_2, mu)

kep = zeros(size(MAT_2,1), 6);
for indm = 1:size(MAT_2,1)

    % incoming leg
    rrIN  = MAT_2(indm,(end-9):(end-9+2));
    vvIN  = MAT_2(indm,(end-6):(end-6+2));

    % evaluate the KEPs
    kep(indm,:) = car2kep([rrIN, vvIN], mu);

end

end