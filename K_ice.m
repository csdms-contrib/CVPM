 function Ki = K_ice(T)

% Finds the thermal conductivity of ice at temperature T based on Yen's (1981)
% expression. Range of validity is 60 K < Tk <= 273.15 K.
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

%   T  = temperature (C)                (vector)
%   Ki = thermal conductivity of ice    (vector)
% ______________________________________________

 Tk = T + 273.15;

% parameters

 a1 = 9.828;
 a2 = 0.0057;

% find the conductivity at temperature T

 Ki    = a1 * exp(-a2 * Tk);
 L     = T > 0;
 Ki(L) = a1 * exp(-a2 * 273.15);

 if any(Tk < 60)
   disp(' ')
   disp('K_ice: temperature is beyond range of validity.')
 end
