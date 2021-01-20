function [eliminate] = find_alreadyVisitedNode(path)

% prendi la ultima riga di path, cioè quella appena aggiunta, e vedi se già
% presente; se si, allora non aggiungerla

rpLast = path(end,14);
ELast  = path(end,15);

drps = path(2:end-1, 14) - rpLast;
dEs  = path(2:end-1, 15) - ELast;

idxsDrps = find(abs(drps) < 0.1);
idsEs    = find(abs(dEs)  < 1e-3);

if ~isempty(idxsDrps) && ~isempty(idsEs) && length(idxsDrps) == length(idsEs) && length(double((idxsDrps == idsEs))) == length(idsEs)
    eliminate = 1; % yes
else
    eliminate = 0; % no
end

end