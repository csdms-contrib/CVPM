 function t = convert_tunits(BC_file,tBC,t_unitsBC,t_units)

% Converts the time units found in the BC_file to agree with the problem's
% declared time units.

%   tBC         BC time vector
%   t_unitsBC   time units associated with tBC (found in BC file)
%   t_units     the problem's time units (found in the namelist file)
%   t           time vector in the problem's time units

% Available time units (t_units):

%   'seconds'
%   'hours'
%   'days'
%   'years'
% ______________________________________________

% Written by:

%   Gary Clow
%   Institute of Arctic and Alpine Research
%   University of Colorado
%   Boulder, Colorado USA
%   Email: gary.clow@colorado.edu
% ______________________________________________

 if ~strcmp(t_units,t_unitsBC)

   switch t_units
   case 'seconds'
     switch t_unitsBC
     case 'days'
       t = 86400 * tBC;
     case 'years'
       t = 86400 * 365 * tBC;
     end

   case 'hours'
     switch t_unitsBC
     case 'days'
       t = 24 * tBC;
     case 'years'
       t = 24 * 365 * tBC;
     end

   case 'days'
     switch t_unitsBC
     case 'years'
       t = 365 * tBC;
     end

   case 'years'
     switch t_unitsBC
     case 'days'
       t = tBC / 365;
     end
   end

   disp(' ')
   disp([BC_file ' time units are different from t_units, I"m converting.'])

 else
   t = tBC;
 end
