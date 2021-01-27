function [tempM_2, tempF_1] = nextLeg(MAT_2, MAT_3, gamma, tolDM, mu, planetsList)

tempM_2 = [];
tempF_1 = [];
for indi = 1:size(MAT_2,1)
    indi;
    % extract info from MAT
    rrIN = MAT_2(indi, end-9:end-9+2);
    vvIN = MAT_2(indi, end-6:end-6+2);
    plIN = MAT_2(indi, end-3);
    tIN  = MAT_2(indi, end-2);
    
    [hmin, hmax] = maxmin_flybyAltitude(plIN);
    
    for indgamma = 1:length(gamma)
        
        gam             = gamma(indgamma);
        [rrOU, vvOU]    = flyby_ciccopl(rrIN, vvIN, plIN, hmin, gam, mu);
        [planetsTarget] = generateTargetBodies(planetsList, rrOU, vvOU);
        
        for indj = 1:length(planetsTarget)
            
            plT = planetsTarget(indj);
            
            [X1, FVAL1] = findPeriAltitude_I(rrIN, vvIN, plIN, tIN, gam, plT, hmin, hmax, planetsList);
            [X2, FVAL2] = findPeriAltitude_II(rrIN, vvIN, plIN, tIN, gam, plT, hmin, hmax, planetsList);
            
            if plT <= 4
                
                [X3, FVAL3] = findPeriAltitude_III(rrIN, vvIN, plIN, tIN, gam, plT, hmin, hmax, planetsList);
                [X4, FVAL4] = findPeriAltitude_IV(rrIN, vvIN, plIN, tIN, gam, plT, hmin, hmax, planetsList);
                
            end

            if FVAL1 < tolDM % && X1 ~= hmin
                [~, rrIN_21, vvIN_21, tPL21] = Intersection_I(X1, gam, rrIN, vvIN, plIN, tIN, plT, planetsList);
                tempM_2 = [tempM_2; [MAT_2(indi,:) [rrIN_21 vvIN_21 plT tPL21 1 0]]];
                tempF_1 = [tempF_1; [MAT_3(indi,:) gam X1]];
                
                [path] = constructPath(tempM_2(end,:), tempF_1(end,:), mu);
                [eliminate] = find_alreadyVisitedNode(path); % elimina i nodi già visitati
                if eliminate == 1
                    tempM_2(end,:) = [];
                    tempF_1(end,:) = [];
                end
            end

            if FVAL2 < tolDM % && X2 ~= hmin
                [~, rrIN_22, vvIN_22, tPL22] = Intersection_II(X2, gam, rrIN, vvIN, plIN, tIN, plT, planetsList);
                tempM_2 = [tempM_2; [MAT_2(indi,:) [rrIN_22 vvIN_22 plT tPL22 2 0]]];
                tempF_1 = [tempF_1; [MAT_3(indi,:) gam X2]];
                
                [path] = constructPath(tempM_2(end,:), tempF_1(end,:), mu);
                [eliminate] = find_alreadyVisitedNode(path); % elimina i nodi già visitati
                if eliminate == 1
                    tempM_2(end,:) = [];
                    tempF_1(end,:) = [];
                end
            end

            if plT <= 4

                if FVAL3 < tolDM % && X3 ~ hmin
                    [~, rrIN_23, vvIN_23, tPL23] = Intersection_III(X3, gam, rrIN, vvIN, plIN, tIN, plT, planetsList);
                    tempM_2 = [tempM_2; [MAT_2(indi,:) [rrIN_23 vvIN_23 plT tPL23 3 0]]];
                    tempF_1 = [tempF_1; [MAT_3(indi,:) gam X3]];
                    
                [path] = constructPath(tempM_2(end,:), tempF_1(end,:), mu);
                [eliminate] = find_alreadyVisitedNode(path); % elimina i nodi già visitati
                if eliminate == 1
                    tempM_2(end,:) = [];
                    tempF_1(end,:) = [];
                end
                end

                if FVAL4 < tolDM % && X4 ~ hmin
                    [~, rrIN_24, vvIN_24, tPL24] = Intersection_IV(X4, gam, rrIN, vvIN, plIN, tIN, plT, planetsList);
                    tempM_2 = [tempM_2; [MAT_2(indi,:) [rrIN_24 vvIN_24 plT tPL24 4 0]]];
                    tempF_1 = [tempF_1; [MAT_3(indi,:) gam X4]];
                    
                [path] = constructPath(tempM_2(end,:), tempF_1(end,:), mu);
                [eliminate] = find_alreadyVisitedNode(path); % elimina i nodi già visitati
                if eliminate == 1
                    tempM_2(end,:) = [];
                    tempF_1(end,:) = [];
                end
                end

            end
        
        end
        
    end
    
end

end