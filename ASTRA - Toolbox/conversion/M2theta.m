function theta = M2theta(M,e)

% compute the true anomaly from mean anomaly with Newton method
%
% INPUT
% - M : mean anomaly (rad)
% - e : eccentricity
%
% OUTPUT : 
% - theta : true anomaly (rad)

E = M; % initial guess for the solution

if e<1 % elliptical orbits
    while abs((M-E + e*sin(E)))>1e-13
        ddf = (e*cos(E)-1);
        E  = E-(M-E + e*sin(E))/ddf;
    end
    
    % from eccentric anomaly to true anomaly
    theta = 2*atan2(sqrt(1+e)*sin(E/2), sqrt(1-e)*cos(E/2));    
    
else % parabolic and hyperbolic transfers 
    for i=1:20
        ddf = (1 - e*cosh(E));
        E  = E-(M + E - e*sinh(E))/ddf;
        
        % from eccentric anomaly to true anomaly
        theta  = 2*atan(sqrt((1+e)/(e-1))*tanh(E/2));
    end
end
end
