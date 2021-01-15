% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  function [planet] = planet_names_GA(planet_id)
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%{
  This function returns the name of the month and the planet
  corresponding, respectively, to the numbers "month_id" and
  "planet_id".
 
  months    - a vector containing the names of the 12 months
  planets   - a vector containing the names of the 9 planets
  month_id  - the month number (1 - 12)
  planet_id - the planet number (1 - 9)
 
  User M-functions required: none
%}
% ------------------------------------------------------------------


planets = ['Y'
           'V'
           'E'
           'M'
           'J'
           'S'
           'U'
           'N'
           'P'];

planet  = planets(planet_id, 1);
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
end 