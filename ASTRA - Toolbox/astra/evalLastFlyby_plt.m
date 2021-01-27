function [MAT_2, MAT_3] = evalLastFlyby_plt(MAT_2, MAT_3, plt, gamma, planetsList)

% evalLastFlyby_plt : this is the version for plt

mu = 132724487690;

inditokeep = [];

for indi = 1:size(MAT_2,1)
   
    % extract info from MAT
    rrIN = MAT_2(indi, end-9:end-9+2);
    vvIN = MAT_2(indi, end-6:end-6+2);
    plIN = MAT_2(indi, end-3);
    
    % evaluate the min. altitude flyby
    [hmin] = maxmin_flybyAltitude(plIN);
    
    for indgamma = 1:length(gamma)
        
        gam = gamma(indgamma);
        [rrOU, vvOU]    = flyby_ciccopl(rrIN, vvIN, plIN, hmin, gam, mu);
        [planetsTarget] = generateTargetBodies(planetsList, rrOU, vvOU);
                
        if ~isempty(find(planetsTarget == plt))
            inditokeep = [inditokeep; indi];
            break
        end
    end
    
end

MAT_2 = MAT_2(inditokeep,:);
MAT_3 = MAT_3(inditokeep,:);

end