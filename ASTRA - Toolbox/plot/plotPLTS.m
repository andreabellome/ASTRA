function [fig] = plotPLTS(pl)

mu = 132724487690;
AU = 149597870.7;

th = linspace(0, 2*pi);

fig = figure;
hold on; grid on; axis equal;
xlabel('x - AU'); ylabel('y - AU');

for indi = 1:length(pl)

    PL = pl(indi);
    
    rr = ephj2000(PL, 0);
    xx = norm(rr).*cos(th);
    yy = norm(rr).*sin(th);
    
    plot(xx./AU, yy./AU, 'k', 'linewidth', 0.5, 'handlevisibility', 'off');
    
end

end