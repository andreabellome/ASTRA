function M = theta2M(theta, e)

% this function computes the mean anomaly from true anomaly and eccentricity
% 
% INPUT : 
% - theta : true anomaly (rad)
% - e     : eccentricity
% 
% OUTPUT:
% - M : mean anomaly (rad)

if e < 1 % elliptical orbits
E = 2*atan(sqrt((1-e)/(1+e))*tan(theta/2)); % from true to eccentric anomaly
M =  E - e*sin(E); % from eccentric to mean anomaly
elseif e > 1 % parabolic and hyperbolic
E = 2*atanh(sqrt((e-1)/(e+1))*tan(theta/2));  % from true to eccentric anomaly
M = e*sinh(E) - E; % from eccentric to mean anomaly
end

end
