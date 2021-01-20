function [DM] = meanAnDifference(M1, M2)

% this function computes the anomaly difference between M1 and M2
%
% INPUT : 
% - M1 : first anomaly (rad)
% - M2 : second anomaly (rad)
% 
% OUTPUT : 
% - DM : anomaly difference (rad)

% TBR (Andrea) -> maybe letting DM to be negative would be beneficial for non linear solver

if abs(M1 - M2) > pi
    DM = 2*pi - abs(M1 - M2);
else
    DM = abs(M1 - M2);
end
    
end
