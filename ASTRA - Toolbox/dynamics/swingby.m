function v2 = swingby(v1, rp, mu, gam, n_r)

% swingby.m - Outgoing velocity after the swingby of a planet.
%
% PROTOTYPE:
%   v2 = swingby(v1, rp, mu, gamma, n_r);
%
% DESCRIPTION:
%   Computes the outgoing velocity after the swingby of a planet, given
%   the incoming velocity and swingby parameters.
%
% INPUT:
%   v1[3]       Incoming velocity, before the swingby, relative to the
%               planet [L/T].
%   rp[1]       Radius of pericentre of the hyperbola [L].
%   mu [1]      Gravitational constant of the planet [L^3/(M*T^2)].
%   gamma[1]	Plane angle [rad]. This angle identifies the inclination of
%               the hyperbola plane around the incoming velocity vector,
%               and is the angle between the vector n_r and the vector
%               normal to the hyperbola plane.
%   n_r[3]      Reference vector, used as a origin to measure gamma. In
%               principle, this vector is arbitrary.
%               A choice can be to use the normal to the plane containing
%               the incoming velocity and the heliocentric velocity of the
%               planet: n_r = cross(v1, vp) / norm(cross(v1, vp));
%               being vp the heliocentric velocity of the planet.
%
% OUTPUT:
%   v2[3,1]     Outgoing velocity, after the swingby, relative to the
%               planet [L/T].
%
% CALLED FUNCTIONS:
%   eulerAxisAngle

mu_rp     = mu/rp;
theta_inf = acos(-mu_rp/(norm(v1)^2+mu_rp)); % See Kaplan "Modern spacecraft dynamics and control", pag. 93
delta     = 2*theta_inf-pi;                  % Deflection angle

n_pi = eulerAxisAngle(n_r, v1, gam);        % Rotates n_r around v1
v2   = eulerAxisAngle(v1, n_pi, delta);     % Rotates v1 around n_pi

return