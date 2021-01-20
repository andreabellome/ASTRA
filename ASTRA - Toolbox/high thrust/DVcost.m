function [dv] = DVcost(kep, kept)

mu = 132724487690;

% given two orbits, find the dv between them
a1 = kep(1); e1 = kep(2);
a2 = kept(1); e2 = kept(2);

rp1 = a1*(1 - e1);
E1  = -mu/(2*a1);

rp2 = a2*(1 - e2);
E2  = -mu/(2*a2);

% find the initial orbit
a1 = -mu/(2*E1);
e1  = 1 - rp1/a1;
ra1 = a1*(1 + e1);

% periapsis/apoapsis velocity of initial orbit
vp1 = sqrt((2*mu)/(ra1 + rp1)*(ra1/rp1));
va1 = sqrt((2*mu)/(ra1 + rp1)*(rp1/ra1));

% find the arrival orbit
a2 = -mu/(2*E2);
e2  = 1 - rp2/a2;
ra2 = a2*(1 + e2);

% periapsis/apoapsis velocity of arrival orbit
vp2 = sqrt((2*mu)/(ra2 + rp2)*(ra2/rp2));
va2 = sqrt((2*mu)/(ra2 + rp2)*(rp2/ra2));

% the ants should not go for hyperbolic transfers
if (e1 >= 1) && (e2 ~= 1)
    dv = 1e99;
elseif (e2 >= 1) && (e1 ~= 1)
    dv = 1e99;
elseif (e1 >= 1) && (e2 >= 1)
    dv = 1e99;
    
else

if (rp1<rp2) && (ra1<ra2)
    
    vpt = sqrt((2*mu)/(rp1 + ra2)*(ra2/rp1));
    vat = sqrt((2*mu)/(rp1 + ra2)*(rp1/ra2));
    
    dv1 = vpt - vp1;
    dv2 = va2 - vat;
    
    dv = dv1 + dv2;
    
elseif (rp1<rp2) && (ra1>ra2)
    
    vpt = sqrt((2*mu)/(rp2 + ra1)*(ra1/rp2));
    vat = sqrt((2*mu)/(rp2 + ra1)*(rp2/ra1));
    
    dv1 = vat - va1;
    dv2 = vpt - vp2;
    
    dv = dv1 + dv2;
    
elseif (rp1==rp2) && (ra1<ra2)
    
    dv = vp2 - vp1;
    
elseif (rp1<rp2) && (ra1==ra2)
    
    dv = va2 - va1;
    
elseif (rp1==rp2) && (ra1>ra2)
    
    dv = vp1 - vp2;
    
elseif (rp1>rp2) && (ra1==ra2)
    
    dv = va1 - va2;
    
elseif (rp1>rp2) && (ra1>ra2)
    
    vpt = sqrt((2*mu)/(rp1 + ra2)*(ra2/rp1));
    vat = sqrt((2*mu)/(rp1 + ra2)*(rp1/ra2));
    
    dv1 = vp1 - vpt;
    dv2 = vat - va2;
    
    dv = dv1 + dv2;
    
elseif (rp1>rp2) && (ra1<ra2)
    
    vpt = sqrt((2*mu)/(rp1 + ra2)*(ra2/rp1));
    vat = sqrt((2*mu)/(rp1 + ra2)*(rp1/ra2));
    
    dv1 = vpt - vp1;
    dv2 = vat - va2;
    
    dv = dv1 + dv2;
    
elseif (rp1==rp2) && (ra1==ra2)
    
    dv = 1e-99;
    
end

end

end