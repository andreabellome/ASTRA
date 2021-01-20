function [rpsc, Esc] = SCorbit(alpha, vInf, vPL, rPL)

mu = 132724487690;

asc  = rPL/(2 - ((vInf/vPL)^2 + 1 + 2*(vInf/vPL)*cos(alpha)));
Esc  = -mu/(2*asc);
vsc  = sqrt(vInf^2 + vPL^2 + 2*vInf*vPL*cos(alpha));
beta = atan2((vInf*sin(alpha)),(vInf*cos(alpha) + vPL));
hsc  = rPL*vsc*cos(beta); 

psc  = hsc^2/mu;
esc  = sqrt(1 - psc/asc);
rpsc = asc*(1 - esc);

end