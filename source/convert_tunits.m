 function t = convert_tunits(BC_file,tBC,t_unitsBC,t_units)

% Converts the time units found in the BC_file to agree with the problem's
% declared time units.
% ______________________________________________

%	Copyright (C) 2018, Gary Clow

%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, version 3 of the License.

%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License v3.0 for more details.

%   You should have received a copy of the GNU General Public License
%   along with this program.  If not, see <http://www.gnu.org/licenses/>.

%	Developer can be contacted by,
%	email at:
%		gary.clow@colorado.edu
%	paper mail at:
%		Institute of Arctic and Alpine Research
%		University of Colorado
%		Campus Box 450
%		Boulder, CO 80309-0450 USA
% ______________________________________________

% Notation:

%   tBC       = BC time vector
%   t_unitsBC = time units associated with tBC (found in BC file)
%   t_units   = the problem's time units (found in the namelist file)
%   t         = time vector in the problem's time units

% Available time units (t_units):

%   'seconds'
%   'hours'
%   'days'
%   'years'
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
