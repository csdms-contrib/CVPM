 function Mw = init_Mw(Mtyp,phi_i,phi_u)

% Initializes the total water mass per unit volume.
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

%   Mtyp  = material type                           (M+1)
%   phi_i = volume fraction of ice                  (M+1)
%   phi_u = volume fraction of unfrozen water       (M+1)
%   Mw    = volumetric water mass                   (M+1)
% ______________________________________________

% parameters

 rhoi = 917;            % density of ice
 rhou = 1000;           % density of unfrozen water

% pre-allocate the array

 Mw = zeros(size(Mtyp));

% reset the water mass for rocks, soils, and organic-rich materials

 L = (Mtyp >= 4) & (Mtyp <= 25);

 if any(L)
   Mw(L) = rhoi * phi_i(L) + rhou * phi_u(L);
 end
