 function rho = init_rho(Mtyp,rhom,phi,phi_i,phi_u)

% Initializes the bulk density.
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

%   Mtyp  = material type                               (M+1)
%   rhom  = density of matrix material                  (M+1)
%   phi   = porosity                                    (M+1)
%   phi_i = volume fraction of ice                      (M+1)
%   phi_u = volume fraction of unfrozen water           (M+1)
%   rho   = bulk density                                (M+1)

% Notes:

%   (1) This is an estimate of rho based on the assumption that all the
%       water is in the ice phase.  When init_rho.m is called, we really
%       don't know how much water is ice and how much is liquid.
% ______________________________________________

% parameters

 rhoi =  917;       % density of ice
 rhou = 1000;       % unfrozen water

% pre-allocate the array to the matrix density

 rho = rhom;

% reset bulk density for rocks, soils, and organic-rich layers

 L = (Mtyp >= 4) & (Mtyp <= 25);

 if any(L)

   rho(L) = (1 - phi(L)).* rhom(L) + phi_i(L)*rhoi + phi_u(L)*rhou;
 end
