 function rho = rho_rock(rhom,phi,phi_i,phi_u)

% Finds the density of porous rocks and soils.
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

%   rhom  = matrix density                          (vector)
%   phi   = porosity                                (vector)
%   phi_i = volume fraction of ice                  (vector)
%   phi_u = volume fraction of unfrozen water       (vector)
%   rho   = bulk density                            (vector)
% ______________________________________________

% component densities

 rhoi =  917;       % ice
 rhou = 1000;       % unfrozen water

% find bulk density

 rho = (1 - phi).*rhom + (phi_i)*rhoi + (phi_u)*rhou;
