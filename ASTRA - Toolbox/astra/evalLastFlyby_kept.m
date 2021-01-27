function [savestates, tempM, tempM_2, tempM_3] = evalLastFlyby_kept(MAT_2, kept, tolDV, gamma)

% evalLastFlyby_kept : this is the version for kept

mu = 132724487690;

tempM   = [];
tempM_3 = [];

for indi = 1:size(MAT_2,1)

    tempM_2 = [];
    
    % extract info. from MAT
    rrIN  = MAT_2(indi, end-9:end-9+2);
    vvIN  = MAT_2(indi, end-6:end-6+2);
    plIN  = MAT_2(indi, end-3);
    tIN   = MAT_2(indi, end-2);
    intIN = MAT_2(indi, end-1);
    
    % evaluate the last flyby
    [hmin, hmax] = maxmin_flybyAltitude(plIN);
    h            = linspace(hmin, hmax);

    for indgamma = 1:length(gamma)
        
        [rrOU, vvOU] = flyby_ciccopl_v2(rrIN, vvIN, plIN, h, gamma(indgamma), mu);
        kepOU        = car2kep_v2([rrOU, vvOU], mu);
        if isempty(kepOU(:,2)<1)
            break
        end
        dv           = DVcost_v2(kepOU(kepOU(:,2) < 1,:), kept);
        
        tempM = [tempM;...
            [rrIN.*ones(size(dv,1),3) vvIN.*ones(size(dv,1),3)...
            plIN.*ones(size(dv,1),1) tIN.*ones(size(dv,1),1) intIN.*ones(size(dv,1),1)...
            indi.*ones(size(dv,1),1) gamma(indgamma).*ones(size(dv,1),1) h(kepOU(:,2) < 1)' rrOU(kepOU(:,2) < 1,:) vvOU(kepOU(:,2) < 1) dv]];
        
        tempM_2 = [tempM_2;...
            [rrIN.*ones(size(dv,1),3) vvIN.*ones(size(dv,1),3)...
            plIN.*ones(size(dv,1),1) tIN.*ones(size(dv,1),1) intIN.*ones(size(dv,1),1)...
            indi.*ones(size(dv,1),1) gamma(indgamma).*ones(size(dv,1),1) h(kepOU(:,2) < 1)' rrOU(kepOU(:,2) < 1,:) vvOU(kepOU(:,2) < 1) dv]];
               
    end
    
    [~, row] = min(tempM_2(:,end));
    tempM_3  = [tempM_3; tempM_2(row,:)];
end

idxs = find(tempM_3(:,end) <= tolDV);

if ~isempty(idxs)
    savestates = tempM_3(idxs,:);
else
    savestates = [];
end

end