 function K = K_rock_porous(T,Km,phi,phi_i,phi_u,planet)

% Finds the bulk thermal conductivity of porous rocks and soils at 
% temperature T.
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

%   T     = temperature (C)                         (vector)
%   Km    = mineral grain conductivity              (vector)
%   phi   = porosity                                (vector)
%   phi_i = volume fraction of ice                  (vector)
%   phi_u = volume fraction of unfrozen water       (vector)
%   K     = bulk thermal conductivity               (vector)
% ______________________________________________

% find the relative volume fractions

 psi_i    = phi_i ./ phi;           % ice
 psi_u    = phi_u ./ phi;           % unfrozen water
 psi_a    = 1 - (psi_i + psi_u);    % air
 L        = psi_a < 0;
 psi_a(L) = 0;                      % make sure psi_a doesn't go slightly negative due to numerics

% find the conductivity of the pores

 Kp = K_pore(T,psi_i,psi_u,psi_a,planet);

% find the bulk conductivity of the mineral-matrix/pore mixture

 K = K_BM2r(Km,Kp,phi);
