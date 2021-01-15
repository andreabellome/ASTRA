function [toftot, TOF_MAN] = tofManoeuvres(path, DV, POS, ALIGN)

% perform the last flyby
% - compute the tof towards the first manoeuvre
% - compute the tof towards the second manoeuvre

% eval the post-flyby orbit
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
    
    TOF_MAN = t;
    
else
    
    % TWO MANOEUVRES - HOHMANN
    % prop. until the position of the first manoeuvre
    if POS(1) == 0
        POS(1) = 2*pi;
    end
    while POS(1) < kep(6)
        POS(1) = POS(1) + 2*pi;
    end
    t = kepEq_t(POS(1), kep(1), kep(2), mu, kep(6), tIN*86400);
    t = t - tIN*86400;

    TOF_MAN(1) = t; % first manoeuvre
    
    % apply the first manoeuvre
    [rrpre, vvpre] = FGKepler_dt(kep, t, mu);
    rrpm  = rrpre;
    vvpm  = vvpre + DV(1).*ALIGN(1).*vvpre./norm(vvpre);
    keppm = car2kep([rrpm vvpm], mu);
    
    TOF_MAN(2) = pi*sqrt(keppm(1)^3/mu);
    
end

toftot = sum(TOF_MAN);

end