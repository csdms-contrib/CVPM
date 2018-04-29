 function Kp = K_pore(T,psi_i,psi_u,psi_a,planet)

% Finds the thermal conductivity of the pores within rocks and soils.
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

%   T     = temperature (C)                             (vector)
%   psi_i = relative volume fraction of ice             (vector)
%   psi_u = relative volume fraction of unfrozen water	(vector)
%   psi_a = relative volume fraction of unfrozen air	(vector)
%   Kp    = conductivity of pores                       (vector)
% ______________________________________________

% pre-allocate arrays

 Kp = NaN*ones(size(psi_i));

% find component conductivities

 Ki = K_ice(T);                 % ice
 Ku = K_water(T);               % unfrozen water
 Ka = K_air(T,planet);          % air

% find the conductivity when the pores are filled with ice

 L     = psi_i == 1;
 Kp(L) = Ki(L);

% find the conductivity when the pores are filled with water

 L     = psi_u == 1;
 Kp(L) = Ku(L);

% find the conductivity when the pores are filled with air

 L     = psi_a == 1;
 Kp(L) = Ka(L);

% find the condutivity when there's more than one phase present

 L = (psi_i ~= 1) & (psi_u ~= 1) & (psi_a ~= 1);
 
 Kp(L) = K_pore_3phase(Ki(L),Ku(L),Ka(L),psi_i(L),psi_u(L),psi_a(L));
