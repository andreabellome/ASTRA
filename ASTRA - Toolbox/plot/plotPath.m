function [fig] = plotPath(path)

mu = 132724487690;
AU = 149597870.7;

PLTS = path(:,7);
TFBS = path(:,8);
DTS  = [0; diff(TFBS)];

[fig] = plotPLTS(PLTS);
hold on;
plot3(0, 0, 0, 'r*', 'markersize', 10, 'handlevisibility', 'off');

% next legs
for NLEG = 2:size(path,1)
    
    RR  = path(NLEG, 1:3);
    VV  = path(NLEG, 4:6);
    
    options = odeset('abstol', 1e-12, 'reltol', 1e-12);
    [tt,yy] = ode113(@(t,x) dyn_2BP(t, x, mu), [0 -DTS(NLEG)*86400], [RR, VV], options);

    hold on;
    plot3(yy(:,1)./AU, yy(:,2)./AU, yy(:,3)./AU, 'b', 'linewidth', 2, 'handlevisibility', 'off');
    plot3(yy(1,1)./AU, yy(1,2)./AU, yy(1,3)./AU, 'ro', 'markersize', 8, 'handlevisibility', 'off');
    
    if NLEG == 2
        hold on;
        plot3(yy(end,1)./AU, yy(end,2)./AU, yy(end,3)./AU, 'ks', 'markersize', 10, 'handlevisibility', 'off');
    end
    
end

end