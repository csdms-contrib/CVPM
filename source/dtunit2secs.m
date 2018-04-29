 function secs_dtunit = dtunit2secs(t_units)

% Finds the number of seconds per dt_unit.
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
