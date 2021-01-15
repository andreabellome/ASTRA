function sd = dyn_2BP(t,s,mu)

% dyn_2BP.m - Equation of the Dynamics of the Two-Body Problem in Cartesian
%   reference frame.
%
% PROTOTYPE:
%   y = dyn_2BP(t, s, mu)
%
% DESCRIPTION:
%   Computes the derivative of the state vector (position and velocity) for
%   the two-body problem dynamics, in Cartesian reference frame.
%
% INPUT:
%   t[1]    Time [T].
%   s[6]    State vector (position, velocity) in Cartesian reference frame
%           [L L/T].
%   mu[1]   Planetary constant [L^3/(M*T^2)].
%           mu = G*mp with  mp = mass of the attracting body
%                           G  = universal gravity constant.
%
% OUTPUT
%   sd[6,1] Derivative of the state vector [L/T L/T^2]
%
% CALLED FUNCTIONS:
%   (none)
%
% EXAMPLE
%   One Earth orbit:
%       [x0,v0] = ephSS_car(3,0);
%       [t,s] = ode45(@dyn_2BP,[0, 365.25*86400],[x0,v0],[],astroConstants(4));
%       plot3(s(:,1),s(:,2),s(:,3));
%
% AUTHOR:
%   Matteo Ceriotti, 2006,, MATLAB, dyn_2BP.m
%
% PREVIOUS VERSION:
%   Matteo Ceriotti, 2006, MATLAB, d2b_eq.m
%       - Header and function name in accordance with guidlines.
%
% CHANGELOG:
%   06/02/2007, REVISION: Camilla Colombo
%   04/01/2010, Camilla Colombo: Header and function name in accordance
%       with guidlines.
%
% -------------------------------------------------------------------------

sd = [       s(4:6)       ;...
      -mu/norm(s(1:3))^3*s(1:3)];
  
return