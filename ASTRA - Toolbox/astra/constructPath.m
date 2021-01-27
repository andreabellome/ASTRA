function [path] = constructPath(MATROW, FLYBYROW, mu)

% this function construct path variable from MATROW (each row of MAT is a planetary sequence)
% 
% INPUT : 
% - MATROW : row with planetary sequence
% - FLYBYROW : row with flyby details for the given sequence
% - mu : gravitational parameter
% 
% OUTPUT : 
% - path : matrix containing all the info for a given planetary sequence
%          - row 1 and cols 1:6     -> initial state at the departure
%          - row 2:end and cols 1:6 -> incoming state at the planet encounter
%          - row 1:end and cols 7   -> planetary sequence
%          - row 1:end and cols 8   -> epochs of planetary encounters
%          - row 2:end and cols 11  -> flyby inclinations
%          - row 2:end and cols 12  -> flyby altitudes

path = reshape(MATROW, 10, [])';

FLYBY = [];
for INDF = 1:2:length(FLYBYROW)
    FLYBY = [FLYBY; FLYBYROW(INDF:INDF+1)]; 
end
FLYBY(end+1,:) = [0 0];

path       = [path FLYBY];

for indj = 1:size(path,1)
   
    rrSC = path(indj,1:3);
    vvSC = path(indj,4:6);
    
    [rrPL] = ephj2000(path(indj,7), path(indj,8));
    [rPL]  = astroConstantsj2000(path(indj,7));
    vvPL   = sqrt(mu/rPL).*[-sin(MeanAn(rrSC)) cos(MeanAn(rrSC)) 0];

    vInf  = norm(vvSC - vvPL);
    kep   = car2kep([rrSC, vvSC], mu);
    rp    = kep(1)*(1 - kep(2));
    E     = -mu/(2*kep(1));
    
    path(indj,13) = vInf;
    path(indj,14) = rp;
    path(indj,15) = E;
    
    [MPL] = MeanAn(rrPL);
    [MSC] = MeanAn(rrSC);
    path(indj,17) = meanAnDifference(MSC, MPL);
    
end

% save the total time of flight
totalTOFyears = sum(diff(path(:,8)))/365.25;
lastcol       = NaN.*ones(size(path,1), 1);
lastcol(1,1)  = totalTOFyears;

path(:,16)    = lastcol;

end
