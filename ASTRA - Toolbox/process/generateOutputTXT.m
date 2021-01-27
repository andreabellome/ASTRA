function [out] = generateOutputTXT(path, filename)

% from ASTRA path to output text

AU   = 149597870.7;
plts = path(:,7);

currentFolder = pwd;
filename      = ['\' filename '_ASTRA_result.txt'];
title         = [currentFolder filename];

% out = fopen('ASTRA_result.txt','w');
out = fopen(title,'w');

fprintf(out,'\n');

fprintf(out, '     _/_/_/     _/_/_/  _/_/_/_/_/  _/_/_/    _/_/_/ \n');
fprintf(out, '   _/    _/   _/           _/     _/    _/  _/    _/ \n');
fprintf(out, '  _/_/_/_/     _/_/       _/     _/_/_/    _/_/_/_/  \n');
fprintf(out, ' _/    _/         _/     _/     _/    _/  _/    _/   \n');       
fprintf(out, '_/    _/    _/_/_/      _/     _/    _/  _/    _/    \n');     

fprintf(out,'\n');

fprintf(out,'               - ASTRA solution - \n');

fprintf(out,'\n');

fprintf(out,'-------------------------------------------------------------- \n');

fprintf(out,'\n');

fprintf(out,['Departing planet            : ' planetIdToName(plts(1)) '\n']);
fprintf(out,'Distance from the Sun       : %.4f AU \n', astroConstantsj2000(plts(1))/AU);

fprintf(out,'\n');

fprintf(out,'-------------------------------------------------------------- \n');

fprintf(out,'\n');

fprintf(out,['Arrival planet              : ' planetIdToName(plts(end)) '\n']);
fprintf(out,'Distance from the Sun       : %.4f AU \n', astroConstantsj2000(plts(end))/AU);
fprintf(out,'Departing infinity velocity : %.4f km/s \n', path(1,9));
fprintf(out,'Arrival infinity velocity   : %.4f km/s \n', path(end,13));
fprintf(out,'Time of flight              : %.4f years \n', path(1,16));
fprintf(out,'Departing C3 (from Earth)   : %.4f km^2/s^2 \n', path(1,9)^2);
fprintf(out,'Total cost                  : %.4f km/s \n', path(1,9)+path(end,13));

fprintf(out,'\n');

fprintf(out,'-------------------------------------------------------------- \n');

fprintf(out,'\n');

fprintf(out, 'MGA Details : \n');

fprintf(out,'\n');

fprintf(out, 'Swing-by sequence : ');
fprintf(out, generateOutputSequenceTXT(path)); fprintf(out,'\n');
fprintf(out, 'Departing date    : ['); fprintf(out, num2str(floor(mjd20002date(path(1,8))), '%d')); fprintf(out, ']'); fprintf(out,'\n');
fprintf(out, 'Arrival date      : ['); fprintf(out, num2str(floor(mjd20002date(path(end,8))), '%d')); fprintf(out, ']'); fprintf(out,'\n');

DTS = diff(path(:,8));
fprintf(out, ['Time of flight per leg : ' num2str(DTS(1)) ' days \n']);
for indi = 2:size(path,1)-1
fprintf(out, ['                         ' num2str(DTS(indi)) ' days \n']);
end

fprintf(out,'\n');

fprintf(out,'-------------------------------------------------------------- \n');

fprintf(out,'\n');

fclose(out);

end