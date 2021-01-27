function [dv] = DVcost_v2(kep, kept)

% kep is a matrix

mu = 132724487690;

% given two orbits, find the dv between them
a1 = kep(:,1);
e1 = kep(:,2);

% find the initial orbit
rp1 = a1.*(1 - e1);
ra1 = a1.*(1 + e1);

a2 = kept(1);
e2 = kept(2);

rp2 = a2*(1 - e2);
E2  = -mu/(2*a2);

% periapsis/apoapsis velocity of initial orbit
vp1 = sqrt((2.*mu)./(ra1 + rp1).*(ra1./rp1));
va1 = sqrt((2.*mu)./(ra1 + rp1).*(rp1./ra1));

% find the arrival orbit
a2 = -mu/(2*E2);
e2  = 1 - rp2/a2;
ra2 = a2*(1 + e2);

% periapsis/apoapsis velocity of arrival orbit
vp2 = sqrt((2*mu)/(ra2 + rp2)*(ra2/rp2));
va2 = sqrt((2*mu)/(ra2 + rp2)*(rp2/ra2));

dv  = zeros(size(kep, 1), 1);
vpt = zeros(size(kep, 1), 1);
vat = zeros(size(kep, 1), 1);
    
% compute dv
dv(e1 >= 1,:) = 1e99;

idxs1          = find((rp1<rp2) & (ra1<ra2));
vpt(idxs1,:)   = sqrt((2*mu)./(rp1(idxs1,:) + ra2).*(ra2./rp1(idxs1,:)));
vat(idxs1,:)   = sqrt((2*mu)./(rp1(idxs1,:) + ra2).*(rp1(idxs1,:)./ra2));
dv1(idxs1,:)   = vpt(idxs1,:) - vp1(idxs1,:);
dv2(idxs1,:)   = va2 - vat(idxs1,:);
dv(idxs1,:)    = dv1(idxs1,:) + dv2(idxs1,:);

idxs2 = find((rp1<rp2) & (ra1>ra2));
vpt(idxs2,:) = sqrt((2*mu)./(rp2 + ra1(idxs2,:)).*(ra1(idxs2,:)./rp2));
vat(idxs2,:) = sqrt((2*mu)./(rp2 + ra1(idxs2,:)).*(rp2./ra1(idxs2,:)));
dv1(idxs2,:) = vat(idxs2,:) - va1(idxs2,:);
dv2(idxs2,:) = vpt(idxs2,:) - vp2;
dv(idxs2,:)  = dv1(idxs2,:) + dv2(idxs2,:);

idxs3 = find((rp1==rp2) & (ra1<ra2));
dv(idxs3,:) = vp2 - vp1(idxs3,:);

idxs4 = find((rp1<rp2) & (ra1==ra2));
dv(idxs4,:) = va2 - va1(idxs4,:);

idxs5 = find((rp1==rp2) & (ra1>ra2));
dv(idxs5,:) = vp1(idxs5,:) - vp2;

idxs6 = find((rp1>rp2) & (ra1==ra2));
dv(idxs6,:) = va1(idxs6,:) - va2;

idxs7 = find((rp1>rp2) & (ra1>ra2));
vpt(idxs7,:) = sqrt((2*mu)./(rp1(idxs7,:) + ra2).*(ra2./rp1(idxs7,:)));
vat(idxs7,:) = sqrt((2*mu)./(rp1(idxs7,:) + ra2).*(rp1(idxs7,:)./ra2));
dv1(idxs7,:) = vp1(idxs7,:) - vpt(idxs7,:);
dv2(idxs7,:) = vat(idxs7,:) - va2;
dv(idxs7,:) = dv1(idxs7,:) + dv2(idxs7,:);

idxs8 = find((rp1>rp2) & (ra1<ra2));
vpt(idxs8,:) = sqrt((2*mu)./(rp1(idxs8,:) + ra2).*(ra2./rp1(idxs8,:)));
vat(idxs8,:) = sqrt((2*mu)./(rp1(idxs8,:) + ra2).*(rp1(idxs8,:)./ra2));
dv1(idxs8,:) = vpt(idxs8,:) - vp1(idxs8,:);
dv2(idxs8,:) = vat(idxs8,:) - va2;
dv(idxs8,:)  = dv1(idxs8,:) + dv2(idxs8,:);

idxs9 = find((rp1==rp2) & (ra1==ra2));
dv(idxs9,:) = 1e-99;
    
end