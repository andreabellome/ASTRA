function [fig] = plotSunDist(path)

global y

mu = 132724487690;
AU = 149597870.7;

TFBS = path(:,8);
DTS  = [0; diff(TFBS)];

fig = figure;
hold on;
xlabel('Time from departure [years]'); ylabel('Spacecraft distance from the Sun [AU]');

% next legs
tend(1) = 0;
for NLEG = 2:size(path,1)-1
    
    RR  = path(NLEG, 1:3);
    VV  = path(NLEG, 4:6);
    
    options  = odeset('abstol', 1e-12, 'reltol', 1e-12);
    [~, yy]  = ode113(@(t,x) dyn_2BP(t, x, mu), [0 -DTS(NLEG)*86400], [RR, VV], options);

    [tt, yy] = ode113(@(t,x) dyn_2BP(t, x, mu), [0 DTS(NLEG)*86400], yy(end,:), options);
    
    r  = zeros(1, size(yy,1));
    rp = zeros(1, size(yy,1));
    ra = zeros(1, size(yy,1));
    for indi = 1:size(yy,1)
        kep     = car2kep(yy(indi,:), mu);
        r(indi) = norm(yy(indi,1:3));
        rp(indi) = kep(1)*(1 - kep(2));
        ra(indi) = kep(1)*(1 + kep(2));
    end
    
    hold on;
    plot((tt+sum(tend))./86400/365.25, r./AU, 'k');
    
    tend(NLEG) = tt(end);
    
end

%%%%%%%
% perform the last flyby and propagate until the perihelion
rrIN  = path(end-1,1:3);
vvIN  = path(end-1,4:6);
plIN  = path(end-1,7);
% tIN   = path(end-1,8);
gamma = path(end-1,11);
hp    = path(end-1,12);

[rrOU, vvOU] = flyby_ciccopl(rrIN, vvIN, plIN, hp, gamma, mu);
kep = car2kep([rrOU, vvOU], mu);

% prop. until the perihelion
t = kepEq_t(2*pi, kep(1), kep(2), mu, kep(6), 0);

[tt, yy] = ode113(@(t,x) dyn_2BP(t, x, mu), [0 t], [rrOU, vvOU], options);
r  = zeros(1, size(yy,1));
rp = zeros(1, size(yy,1));
ra = zeros(1, size(yy,1));
for indi = 1:size(yy,1)
    kep      = car2kep(yy(indi,:), mu);
    r(indi)  = norm(yy(indi,1:3));
    rp(indi) = kep(1)*(1 - kep(2));
    ra(indi) = kep(1)*(1 + kep(2));
end

hold on;
plot((tt+sum(tend))./86400/365.25, r./AU, 'k');
%%%%%%%

%%%%%%%
% include vertical lines for flybys
y = get(gca,'ylim');
ylim(y);
for NLEG = 2:size(path,1)-1
    plGA = planet_names_GA(path(NLEG,7));
    vline((TFBS(NLEG)-path(1,8))/365.25, 'k', plGA);
end

end