function [fig] = plotPathManoeuvres(path, POS, kept)

% plot the traj. until the last flyby
[fig] = plotPath(path);

% eval the post-flyby orbit
AU = 149597870.7;
mu = 132724487690;

rrIN  = path(end,1:3);
vvIN  = path(end,4:6);
plIN  = path(end,7);
gamma = path(end,11);
hp    = path(end,12);

tIN   = path(end,8);

[rrOU, vvOU] = flyby_ciccopl(rrIN, vvIN, plIN, hp, gamma, mu);
kep          = car2kep([rrOU, vvOU], mu);

if length(POS) == 1
    % ONLY ONE MANOEUVRE - HOHMANN
    if POS(1) == 0
        POS(1) = 2*pi;
    end
    while POS(1) < kep(6)
        POS(1) = POS(1) + 2*pi;
    end
    % prop. until the position of the manoeuvre
    t = kepEq_t(POS(1), kep(1), kep(2), mu, kep(6), tIN*86400);
    t = t - tIN*86400;

    % plot the prop. traj. until the position of the manoeuvre
    options = odeset('abstol', 1e-12, 'reltol', 1e-12);
    [~, yy] = ode113(@(t,x) dyn_2BP(t, x, mu), [0 t], [rrOU, vvOU], options);

    hold on;
    plot(yy(:,1)./AU, yy(:,2)./AU, 'b', 'linewidth', 2, 'handlevisibility', 'off');
    plot(yy(end,1)./AU, yy(end,2)./AU, 'x', 'markersize', 10, 'linewidth', 3);

    [~, ~, ~, DV, ALIGN] = DVcost_v3(kep, kept);
    
    rrpm  = yy(end,1:3);
    vvpm  = yy(end,4:6) + DV.*ALIGN.*yy(end,4:6)./norm(yy(end,4:6));
    keppm = car2kep([rrpm vvpm], mu);
    tofpm = 2*pi*sqrt(keppm(1)^3/mu);
    
    % plot the target orbit (one period)
    [~, yy] = ode113(@(t,x) dyn_2BP(t, x, mu), [0 tofpm], [rrpm vvpm], options);
    
    hold on;
    plot(yy(:,1)./AU, yy(:,2)./AU, 'm', 'linewidth', 2);
    
    legend('Manoeuvre location', 'Final orbit');
    legend('location', 'northwest');

else
    % TWO MANOEUVRES - HOHMANN
    
    % prop. until the position of the first manoeuvre
    if POS(1) == 0
        POS(1) = 2*pi;
    end
    t = kepEq_t(POS(1), kep(1), kep(2), mu, kep(6), tIN*86400);
    t = t - tIN*86400;

    % plot the prop. traj. until the position of the manoeuvre
    options = odeset('abstol', 1e-12, 'reltol', 1e-12);
    [~, yy] = ode113(@(t,x) dyn_2BP(t, x, mu), [0 t], [rrOU, vvOU], options);

    hold on;
    plot(yy(:,1)./AU, yy(:,2)./AU, 'b', 'linewidth', 2, 'handlevisibility', 'off');
    plot(yy(end,1)./AU, yy(end,2)./AU, 'x', 'markersize', 10, 'linewidth', 3);

    [~, ~, ~, DV, ALIGN] = DVcost_v3(kep, kept);
    % apply the first manoeuvre
    rr     = yy(end,1:3);
    vvPost = yy(end,4:6) + DV(1).*ALIGN(1).*yy(end,4:6)./norm(yy(end,4:6));
    
    kepPost = car2kep([rr, vvPost], mu);
    
    tof_t   = pi*sqrt(kepPost(1)^3/mu);
    [~, yy] = ode113(@(t,x) dyn_2BP(t, x, mu), [0 tof_t], [rr, vvPost], options);

    plot(yy(:,1)./AU, yy(:,2)./AU, 'b', 'linewidth', 2, 'handlevisibility', 'off');
    plot(yy(end,1)./AU, yy(end,2)./AU, 'x', 'markersize', 10, 'linewidth', 3);
    
    % apply the second manoeuvre
    rr     = yy(end,1:3);
    vvPost = yy(end,4:6) + DV(2).*ALIGN(2).*yy(end,4:6)./norm(yy(end,4:6));
    
    kepPost = car2kep([rr, vvPost], mu);
    
    tof_t   = pi*sqrt(kepPost(1)^3/mu);
    [~, yy] = ode113(@(t,x) dyn_2BP(t, x, mu), [0 tof_t], [rr, vvPost], options);

    plot(yy(:,1)./AU, yy(:,2)./AU, 'b', 'linewidth', 2, 'handlevisibility', 'off');

    % plot the final orbit
    [~, yy] = ode113(@(t,x) dyn_2BP(t, x, mu), [0 2*tof_t], [rr, vvPost], options);
    plot(yy(:,1)./AU, yy(:,2)./AU, 'm', 'linewidth', 2);
    
    legend('First manoeuvre', 'Second manoeuvre', 'Final orbit');
    legend('location', 'northwest');
    
end

end