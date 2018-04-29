 function rho = update_rho(Mtyp,rhom,phi,phi_i,phi_u)

% Updates the bulk density field.
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

%   Mtyp  = material type
%   rhom  = matrix density
%   phi   = porosity
%   phi_i = volume fraction of ice
%   phi_u = volume fraction of unfrozen water

% Notes:
%   (1) We only need to do this for porous materials.
% ______________________________________________

% default

 rho = rhom;

% reset bulk density of layers consisting of rocks, soils, or organic-rich layers

 L = (Mtyp >= 4) & (Mtyp <= 25);

 rho(L) = rho_rock(rhom(L),phi(L),phi_i(L),phi_u(L));
