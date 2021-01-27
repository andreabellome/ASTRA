function [pos, vel, posICRF, velICRF] = EphAA_car(pl, t)

% TITLE : Keplerian Elements for Approximate Positions of the Major Planets
% DESCRIPTION :
% Ephemerides with respect to mean ecliptic and equinox of J2000.0, valid
% for the time interval 1800 AD - 2050 AD
%
% INPUT :
% pl = [1 2 3 4 5 6 7 8 9]
% t is MJD2000
%
% OUTPUT :
% pos     : (1x3) position vector of the planet at t
% vel     : (1x3) velocity vector of the planet at t
% posICRF : (1x3) position vector of the planet at t (ICRF frame)
% velICRF : (1x3) velocity vector of the planet at t

% constants of motion
mu = 132724487690;
AU = 149597870.7;

% compute centuries past J2000.0
T = t/36525;

% obliquity at J2000
obl = deg2rad(23.43928); 

% keplerian elements and their rates with respect to mean ecliptic and
% equinox of J2000.0, valid for the time interval 1800 AD - 2050 AD
keppl = [                   0.38709927,      0.20563593,      7.00497902,      252.25032350,     77.45779628,     48.33076593;
                            0.72333566,      0.00677672,      3.39467605,      181.97909950,    131.60246718,     76.67984255;
                            1.00000261,      0.01671123,     -0.00001531,      100.46457166,    102.93768193,      0.0; 
                            1.52371034,      0.09339410,      1.84969142,       -4.55343205,    -23.94362959,     49.55953891; 
                            5.20288700,      0.04838624,      1.30439695,       34.39644051,     14.72847983,    100.47390909; 
                            9.53667594,      0.05386179,      2.48599187,       49.95424423,     92.59887831,    113.66242448;
                            19.18916464,      0.04725744,      0.77263783,      313.23810451,    170.95427630,     74.01692503;
                            30.06992276,      0.00859048,      1.77004347,      -55.12002969,     44.96476227,    131.78422574;
                            39.48211675,      0.24882730,     17.14001206,      238.92903833,    224.06891629,    110.30393684 ];


kepdotpl = [             0.00000037,      0.00001906,     -0.00594749,   149472.67411175,      0.16047689,     -0.12534081;  
                         0.00000390,     -0.00004107,     -0.00078890,    58517.81538729,      0.00268329,     -0.27769418;  
                         0.00000562,     -0.00004392,     -0.01294668,    35999.37244981,      0.32327364,      0.0;           
                         0.00001847,      0.00007882,     -0.00813131,    19140.30268499,      0.44441088,     -0.29257343;  
                        -0.00011607,     -0.00013253,     -0.00183714,     3034.74612775,      0.21252668,      0.20469106; 
                        -0.00125060,     -0.00050991,      0.00193609,     1222.49362201,     -0.41897216,     -0.28867794; 
                        -0.00196176,     -0.00004397,     -0.00242939,      428.48202785,      0.40805281,      0.04240589; 
                         0.00026291,      0.00005105,      0.00035372,      218.45945325,     -0.32241464,     -0.00508664;  
                        -0.00031596,      0.00005170,      0.00004818,      145.20780515,     -0.04062942,     -0.01183482 ];

% STEP 1 : compute the value of each of that planet's six elements
KEP  = keppl(pl,:);
KEPd = kepdotpl(pl,:);

KEPt(1) = KEP(1) + T*KEPd(1); % (au) semi_major_axis
KEPt(2) = KEP(2) + T*KEPd(2); % ( ) eccentricity
KEPt(3) = KEP(3) + T*KEPd(3); % (°) inclination
KEPt(4) = KEP(4) + T*KEPd(4); % (°) mean_longitude
KEPt(5) = KEP(5) + T*KEPd(5); % (°) longitude_of_periapsis
KEPt(6) = KEP(6) + T*KEPd(6); % (°) longitude_of_the_ascending_node

% STEP 2 : compute the argument of perihelion, ω, and the mean anomaly, M
om = KEPt(5) - KEPt(6);
M  = KEPt(4) - KEPt(5);

% STEP 3a : modulus the mean anomaly so that -180° ≤ M ≤ +180°
while M > 180
    M = M - 360;
end

% STEP 3b : obtain the eccentric anomaly, E, from the solution of Kepler's equation
E = M + KEPt(2)*180/pi*sin(M*pi/180); % initial guess
for indi = 1:6
   E = E + (M - (E - KEPt(2)*180/pi*sin(E*pi/180)))/(1 - KEPt(2)*cos(E*pi/180)); 
end

% STEP 4 : compute the planet's heliocentric coordinates in its orbital plane, r', with the x'-axis aligned from the focus to the perihelion
om  = om*pi/180;
OM  = KEPt(6)*pi/180;
E   = E*pi/180;
th  = 2*atan2(sqrt(1+KEPt(2))*sin(E/2), sqrt(1-KEPt(2))*cos(E/2));
i   = KEPt(3)*pi/180;
x0  = KEPt(1)*(cos(E) - KEPt(2));
y0  = KEPt(1)*sqrt(1-KEPt(2)*KEPt(2))*sin(E);
vx0 = sqrt(mu/(KEPt(1)*AU*(1 - KEPt(2)^2)))*(KEPt(2)*sin(th)*cos(th) - (1 + KEPt(2)*cos(th))*sin(th));
vy0 = sqrt(mu/(KEPt(1)*AU*(1 - KEPt(2)^2)))*(KEPt(2)*sin(th)*sin(th) + (1 + KEPt(2)*cos(th))*cos(th));

% STEP 5a : compute the coordinates in the J2000 ecliptic plane, with the x-axis aligned toward the equinox
pos    = zeros(1, 3);
vel    = zeros(1, 3);
pos(1) = (cos(om)*cos(OM) - sin(om)*sin(OM)*cos(i))*x0 + (-sin(om)*cos(OM) - cos(om)*sin(OM)*cos(i))*y0;
pos(2) = (cos(om)*sin(OM) + sin(om)*cos(OM)*cos(i))*x0 + (-sin(om)*sin(OM) + cos(om)*cos(OM)*cos(i))*y0;
pos(3) = (sin(om)*sin(i))*x0 + (cos(om)*sin(i))*y0;
pos    = pos.*AU;
vel(1) = (cos(om)*cos(OM) - sin(om)*sin(OM)*cos(i))*vx0 + (-sin(om)*cos(OM) - cos(om)*sin(OM)*cos(i))*vy0;
vel(2) = (cos(om)*sin(OM) + sin(om)*cos(OM)*cos(i))*vx0 + (-sin(om)*sin(OM) + cos(om)*cos(OM)*cos(i))*vy0;
vel(3) = (sin(om)*sin(i))*vx0 + (cos(om)*sin(i))*vy0;

% STEP 5b : compute the coordinates in the ICRF plane, inclined by obl.
% around x
posICRF = zeros(1,3);
posICRF(1) = pos(1);
posICRF(2) = cos(obl)*pos(2) - sin(obl)*pos(3);
posICRF(3) = sin(obl)*pos(2) + cos(obl)*pos(3);

velICRF    = zeros(1,3);
velICRF(1) = vel(1);
velICRF(2) = cos(obl)*vel(2) - sin(obl)*vel(3);
velICRF(3) = sin(obl)*vel(2) + cos(obl)*vel(3);

end