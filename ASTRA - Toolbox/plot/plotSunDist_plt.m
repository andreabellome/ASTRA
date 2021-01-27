function [fig] = plotSunDist_plt(path)

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
options = odeset('abstol', 1e-12, 'reltol', 1e-12);
for NLEG = 2:size(path,1)-1
    
    RR  = path(NLEG, 1:3);
    VV  = path(NLEG, 4:6);
    
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

% propagate the last leg
rrINend = path(end,1:3);
vvINend = path(end,4:6);
DTSend  = DTS(end);

[~, yy]  = ode113(@(t,x) dyn_2BP(t, x, mu), [0 -DTSend*86400], [rrINend, vvINend], options);
[tt, yy] = ode113(@(t,x) dyn_2BP(t, x, mu), [0 DTSend*86400], yy(end,:), options);
tt = sum(tend) + tt;
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
plot((tt)./86400/365.25, r./AU, 'k');

%%%%%%%
% include vertical lines for flybys
y = get(gca,'ylim');
ylim(y);
for NLEG = 2:size(path,1)
    plGA = planet_names_GA(path(NLEG,7));
    vline((TFBS(NLEG)-path(1,8))/365.25, 'k', plGA);
end

end