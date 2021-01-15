function M = theta2M(theta, e)

if e < 1
E = 2*atan(sqrt((1-e)/(1+e))*tan(theta/2));
M =  E - e*sin(E);
elseif e > 1
E = 2*atanh(sqrt((e-1)/(e+1))*tan(theta/2));
M = e*sinh(E) - E;
end

end