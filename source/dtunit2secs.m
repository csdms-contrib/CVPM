 function secs_dtunit = dtunit2secs(t_units)

% Finds the number of seconds per dt_unit.
% ______________________________________________

% Written by:

%   Gary Clow
%   Institute of Arctic and Alpine Research
%   University of Colorado
%   Boulder, Colorado USA
%   Email: gary.clow@colorado.edu
% ______________________________________________

 switch t_units
 case 'seconds'
   secs_dtunit = 1;
 case 'hours'
   secs_dtunit = 3600;
 case 'days'
   secs_dtunit = 86400;
 case 'weeks'
   secs_dtunit = 7*86400;
 case 'months'
   secs_dtunit = 365*86400/12;
 case 'years'
   secs_dtunit = 365*86400;
 end
