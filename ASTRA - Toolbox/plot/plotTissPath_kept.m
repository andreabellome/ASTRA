function [fig] = plotTissPath_kept(path)

% perform the flyby at the last planet of the sequence

mu = 132724487690;
AU = 149597870.7;

% plot the Tisserand route;

PLTS = path(:,7);
VINF = path(:,13);
RPS  = path(:,14);
ES   = path(:,15);

fig = figure;
hold on; grid on;
xlabel('r_p - AU'); ylabel('E - km^2/s^2');

% plot flybys
for indj = 1:size(path,1)-1
       
   if indj == 1
       
       hold on;
       plot(RPS(indj)./AU, ES(indj), 'ks', 'markersize', 10, 'handlevisibility', 'off');
       
   elseif indj == 2
       
   else
       
       hold on;
       plot(RPS(indj)./AU, ES(indj), 'k.', 'markersize', 15, 'handlevisibility', 'off');
       
   end
      
end

% plot contours
for indj = 1:size(path,1)
    
   [rpscCONT, EscCONT] = generateContours_CircCopl(PLTS(indj), VINF(indj));
   hold on;
   
   if PLTS(indj) == 2
       color = ['g'];
   elseif PLTS(indj) == 3
       color = ['k'];
   elseif PLTS(indj) == 4
       color = ['m'];
   end
   plot(rpscCONT./AU, EscCONT, color, 'linewidth', 0.5, 'handlevisibility', 'off');
end

hold on;
plot(RPS(end)./AU, ES(end), 'k.', 'markersize', 15, 'handlevisibility', 'off');

if norm(double((PLTS) == 5)) ~= 0
    % metti limiti solo se c'è Giove
    xlim([0 2]);
end

% plot tisserand trajectory
for indi = 2:size(path,1)-1
    
    rr = path(indi,1:3);
    vv = path(indi,4:6);
    
    kep = car2kep([rr vv], mu);
    
    r_p_1 = kep(1)*(1 - kep(2));
    r_a_1 = kep(1)*(1 + kep(2));
    
    % compute in front/behind the planet
    Eprev   = -mu/(2*kep(1));
    nextRow = indi + 1;
    rrnext  = path(nextRow, 1:3);
    vvnext  = path(nextRow, 4:6);
    kepnext = car2kep([rrnext vvnext], mu);
    Enext   = -mu/(2*kepnext(1));
    if Enext > Eprev
        tj = 1;
    else
        tj = 0;
    end
        
    Fb_planet = path(indi,7);
    
    [r_p_pl , r_a_pl ] = Tisserand_Trajectory ( r_p_1 , r_a_1 , tj , Fb_planet , path(indi,12), 100 );
    a_pl = 0.5.*(r_p_pl + r_a_pl);
    E_pl = -mu./(2*a_pl);
    
    hold on;
    plot(r_p_pl(1,:)'./AU , E_pl , 'b' , 'LineWidth' , 2 , 'handlevisibility', 'off');
    
end

rrIN   = path(end,1:3);
vvIN   = path(end,4:6);
plIN   = path(end,7);
gamma  = path(end,11);
hp     = path(end,12);
[rrOU, vvOU] = flyby_ciccopl(rrIN, vvIN, plIN, hp, gamma, mu);
kep          = car2kep([rrOU, vvOU], mu);
rp = kep(1)*(1 - kep(2));
E  = -mu/(2*kep(1));

hold on;
plot(rp./AU, E, 'k.', 'markersize', 15, 'handlevisibility', 'off');

kepIN = car2kep([rrIN, vvIN], mu);
rpIN  = kepIN(1)*(1 - kepIN(2));
EIN   = -mu/(2*kepIN(1));
raIN  = kepIN(1)*(1 + kepIN(2));
if E > EIN
    tj = 1;
else
    tj = 0;
end

[r_p_pl , r_a_pl ] = Tisserand_Trajectory ( rpIN , raIN , tj , plIN , hp, 100 );
a_pl = 0.5.*(r_p_pl + r_a_pl);
E_pl = -mu./(2*a_pl);

hold on;
plot(r_p_pl(1,:)'./AU , E_pl , 'b' , 'LineWidth' , 2 , 'handlevisibility', 'off');

end