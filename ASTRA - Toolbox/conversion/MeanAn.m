function [MPL] = MeanAn(rrPL)

% this function computes the mean anomaly from the position vector
%
% INPUT : 
% - rrPL : planet position vector (km)
%
% OUTPUT : 
% - MPL : mean anomaly (rad)

% read the position vector components
xPL = rrPL(1);
yPL = rrPL(2);

% find mean anomaly and wrap to 2pi
MPL = atan2(yPL, xPL);
MPL = wrapTo2Pi(MPL);

end
