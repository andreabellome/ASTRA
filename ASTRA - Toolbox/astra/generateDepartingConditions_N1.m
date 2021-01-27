function [MMAT_2, DV0, meshTDEP, meshTOF] = ...
    generateDepartingConditions_N1(s, DEP_OPTIONS, ARR_OPTIONS, vinfMin, vinfMax)

% this function generates the first leg for ASTRA toolbox
% (only direct transfers - one revolution on the first leg)

mu = 132724487690;

t0       = DEP_OPTIONS(1); % mission start date
Dt0      = DEP_OPTIONS(2); % mission start date range
tsteptt0 = DEP_OPTIONS(3); % start interval

minTOF   = ARR_OPTIONS(1); % minimum transfer duration range 
maxTOF   = ARR_OPTIONS(2); % maximum transfer duration range 
tstepTT  = ARR_OPTIONS(3); % duration interval

% generate departure and arrival grids
tt0      = t0:tsteptt0:(t0+Dt0);
TT       = minTOF:tstepTT:maxTOF;

DV0             = zeros(length(tt0), length(TT));
meshTDEP        = zeros(length(tt0), length(TT));
meshTOF         = zeros(length(tt0), length(TT));
meshDateMJD2000 = zeros(length(tt0), length(TT));

dep_cond = [];

% initialise the waitbar
f = waitbar(0, '', 'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');

for indi = 1:length(tt0)
    
    TDEP    = tt0(indi);
    
    [r1,v1] = ephj2000(s(1), TDEP);
            
    for indj = 1:length(TT)
        
        TOF     = TT(indj);
        
        % Build meshes for contour plot axis
        meshTDEP(indi,indj)        = TDEP;
        meshTOF(indi,indj)         = TOF;
        meshDateMJD2000(indi,indj) = TDEP + TOF; % this is the arrival date at the planet
        
        [r2] = ephj2000(s(2),TDEP + TOF);
                
        % one revolution
        [~,~,~,~,v1Short, v2Short,~,~] = lambertMR(r1, r2, TOF*86400, mu, 0, 1, 0, 0);
        dvShort = norm(v1Short - v1);
        
        [~,~,~,~,v1Long, v2Long,~,~] = lambertMR(r1, r2, TOF*86400, mu, 0, 1, 1, 0);
        dvLong = norm(v1Long - v1);
        
        [DV0(indi,indj), b] = min([dvShort dvLong]);
        
        % save pairs of states 
        if b == 1
            % short transfer is saved
            rr     = r1;
            vv     = v1Short;
            vInfD  = norm(v1Short - v1);
            alphaD = acos(dot(v1Short - v1,v1)/(vInfD*norm(v1)));
            
            if vInfD < vinfMax
                rrIN     = r2;
                vvIN     = v2Short;
                dep_cond = [dep_cond; [rr vv s(1) TDEP vInfD alphaD rrIN vvIN s(2) (TDEP + TOF) b 0]];
                
            end
           
        elseif b == 2
            % long transfer is saved
            rr     = r1;
            vv     = v1Long;
            vInfD  = norm(v1Long - v1);
            alphaD = acos(dot(v1Long - v1,v1)/(vInfD*norm(v1)));
            
            if vInfD < vinfMax
                rrIN     = r2;
                vvIN     = v2Long;
                dep_cond = [dep_cond; [rr vv s(1) TDEP vInfD alphaD rrIN vvIN s(2) (TDEP + TOF) b 0]];
                
            end
            
        end
        
    end
    
    % update the waitbar
    w = waitbar(indi/length(tt0), f, sprintf('Computing first leg - N1...'));

end
delete(w); % close the waitbar

MAT_2  = dep_cond;
VINF   = MAT_2(:,9);
MAT_2  = MAT_2(find(VINF <= vinfMax & VINF >= vinfMin),:);
MMAT_2 = MAT_2;

end