function [fig] = plotPorkChop(s, DV0, meshTDEP, meshTOF)

minDV0       = min(min(DV0));
V_inf_levels = minDV0:1:50;

N = zeros(size(meshTDEP, 1), size(meshTDEP, 2));
for indi = 1:size(meshTDEP, 1)
    for indj = 1:size(meshTDEP, 2)
        date = mjd20002date(meshTDEP(indi,indj));
        date = date(1:3);
        N(indi,indj) = datenum(date);
    end
end

fig = figure;
col1 = [0.8,0.2,0.2];

hold on;
[c1,h1] = contourf(meshTOF, N, DV0, V_inf_levels, 'color', col1, 'linewidth', 0.5, 'handlevisibility', 'off');
% colorbar;
hcb = colorbar;
title(hcb,'v_{\infty,dep}')
set(gca, 'YDir','reverse');
xlabel('Time of Flight [days]'); ylabel('Departing Date');
if s == [3 2]
    title('Earth - Venus transfer');
else
    title('Earth - Mars transfer');
end

datetick('y',26);

hold off;

end