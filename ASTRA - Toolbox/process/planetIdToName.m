function [name] = planetIdToName(idPL)

names = {'Mercury', 'Venus', 'Earth', 'Mars', 'Jupiter', 'Saturn', 'Uranus', 'Neptune'};
name  = char(names(idPL));

end