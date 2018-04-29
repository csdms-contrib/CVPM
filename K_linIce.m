 function Ki = K_linIce(T,K0)

% Finds the thermal conductivity of a substance for which K is linearly
% dependent on temperature.  For this we choose properties similar to ice.
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

%   T  = temperature (C)                      (vector)
%   K0 = conductivity @ 0 C                   (vector)
%   Ki = conductivity at temperature T        (vector)
% ______________________________________________

% parameters

 T0   = 0;              % reference temperature (C)
 dKdT = -0.0125;        % dK/dT

% find the conductivity

 Ki = K0 + dKdT *(T - T0);
