function [MAT_2, MAT_3, MMMAT_2dv] = apply_beamSearch(bw, MAT_2, MAT_3, kept)

% check the DV from target
dv   = find_dvToTargetFromMMAT(MAT_2, kept);

% apply BS
MMMAT_2dv = [MAT_2 dv];
MMMAT_3dv = [MAT_3 dv];
MMMAT_2dv = sortrows(MMMAT_2dv, size(MMMAT_2dv,2));
MMMAT_3dv = sortrows(MMMAT_3dv, size(MMMAT_3dv,2));

if size(MMMAT_2dv, 1) <= bw

else
    MMMAT_2dv = MMMAT_2dv(1:bw,:);
    MMMAT_3dv = MMMAT_3dv(1:bw,:);

    MAT_2    = MMMAT_2dv(:,1:end-1);
    MAT_3    = MMMAT_3dv(:,1:end-1);
end

end