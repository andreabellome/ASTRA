function [rrf, vvf] = FGKepler_dt(kep1, dt, mu)

% this function computes the final position and velocity vectors given intial keplerian elements,
% transfer time and the gravitational parameter of the central body
%
% INPUT : 
% - kep1 : initial keplerian elements [1x6]
% - dt   : transfer time (secs)
% - mu   : gravitational parameter of the central body (km3/s2)
%
% OUTPUT : 
% - rrf : final position vector (km)
% - vvf : final velocity vector (km)

T    = 2*pi*sqrt(kep1(1)^3/mu); % period of the transfer orbit
n    = 2*pi/T; % mean motion of the transfer orbit
M1   = theta2M(kep1(6), kep1(2)); % initial mean anomaly
M2   = M1 + n*dt; % final mean anomaly
th2  = M2theta(M2,kep1(2)); % from final mean to final true anomaly

kep2 = [kep1(1) kep1(2) kep1(3) kep1(4) kep1(5) th2]; % keplerian elements of the final orbit
car2 = kep2car(kep2,mu); % from keplerian to cartesian elements

rrf = car2(1:3)'; % extract final position vector
vvf = car2(4:6)'; % extract final velocity vector
rrf = rrf'; vvf = vvf'; % make row vectors

end
