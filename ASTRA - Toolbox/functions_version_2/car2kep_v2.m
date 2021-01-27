function [kep] = car2kep_v2(in,mu)

% car2kep
% in can also be a [n,6] matrix

elimit = 0.00000001;

r  = in(:,1:3);
v  = in(:,4:6);
nr = sqrt(r(:,1).^2 + r(:,2).^2 + r(:,3).^2);    % Norm of r

% Angular momentum vector: h = cross(r,v)
h(:,1) = r(:,2).*v(:,3) - r(:,3).*v(:,2);
h(:,2) = r(:,3).*v(:,1) - r(:,1).*v(:,3);
h(:,3) = r(:,1).*v(:,2) - r(:,2).*v(:,1);
nh     = sqrt(h(:,1).^2 + h(:,2).^2 + h(:,3).^2);    % Norm of h

% Inclination
i = acos(h(:,3)./nh);

idxs1 = find(i~=0);
idxs2 = find(i~=pi);
n = [1, 0, 0].*ones(size(h,1),3);
if ~(isempty(idxs1)) && ~(isempty(idxs2))
    n(idxs1,1) = [-h(idxs1,2), h(idxs1,1), zeros(size(h(idxs1,:),1),1)]./sqrt(h(idxs1,1).^2 + h(idxs1,2).^2);
    n(idxs2,1) = [-h(idxs2,2), h(idxs2,1), zeros(size(h(idxs2,:),1),1)]./sqrt(h(idxs2,1).^2 + h(idxs2,2).^2);
end

% Argument of the ascending node
Om          = acos(n(:,1));
idxsn       = find(n(:,2) < 0);
Om(idxsn,:) = mod(2*pi - Om(idxsn,:), 2*pi);

% Parameter
p = nh.^2./mu;

ev(:,1) = 1/mu.*(v(:,2).*h(:,3) - v(:,3).*h(:,2)) - r(:,1)./nr;
ev(:,2) = 1/mu.*(v(:,3).*h(:,1) - v(:,1).*h(:,3)) - r(:,2)./nr;
ev(:,3) = 1/mu.*(v(:,1).*h(:,2) - v(:,2).*h(:,1)) - r(:,3)./nr;
e       = sqrt(ev(:,1).^2 + ev(:,2).^2 + ev(:,3).^2);    % Eccentricity (norm of eccentricity vector)
ne = e;
idxse = find(e < elimit);
ev(idxse,:) = n(idxse,:);
ne(idxse,1) = 1;

om = acos(min(max((n(:,1).*ev(:,1) + n(:,2).*ev(:,2) + n(:,3).*ev(:,3))./ne, -1), 1)); % acos(dot(n,ev)/ne)
idxsom = find(dot(h,cross(n,ev)) < 0 );
om(idxsom,:) = mod(2*pi-om(idxsom,:),2*pi);

% Semi-major axis
a = p./(1-e.^2);

% True anomaly: acos(dot(ev,r)/ne/nr);
th = acos(min(max((ev(:,1).*r(:,1) + ev(:,2).*r(:,2) + ev(:,3).*r(:,3))./ne./nr,-1),1));

dothcrossevr  = h(:,1).*(ev(:,2).*r(:,3) - ev(:,3).*r(:,2))...
    + h(:,2).*(ev(:,3).*r(:,1) - ev(:,1).*r(:,3))...
    + h(:,3).*(ev(:,1).*r(:,2) - ev(:,2).*r(:,1));

idxsth = find(dothcrossevr < 0);
th(idxsth,:) = mod(2*pi-th(idxsth,:),2*pi);

kep = [a, e, i, Om, om, th];

end