function [path] = constructPath(MATROW, FLYBYROW, mu)

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