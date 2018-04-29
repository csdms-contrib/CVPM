 function [rhocp,C] = C_veg(T,Mtyp,rhom,cpm0,phi,phi_i,phi_u,dphiudT)

% Finds the volumetric heat capacity of organic-rich material at temperature T.
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

%   T       = temperature (C)                           (vector)
%   Mtyp    = material type                             (vector)
%   rhom    = mean matrix density                       (vector)
%   cpm0    = pecific heat of matrix at 20 C            (vector)
%   phi     = porosity                                  (vector)
%   phi_i   = volume fraction of ice                    (vector)
%   phi_u   = volume fraction of unfrozen water         (vector)
%   dphiudT = d(phi_u)/dT                               (vector)
%   rhocp   = lattice component of the heat capacity	(vector)
%   C       = volumetric heat capacity                  (vector)
% ______________________________________________

% parameters

 rhow = 1000;           % density of water
 Lf   = 3.34e05;        % latent heat of fusion (water)

% find component heat capacities

 rhocpm = C_organic(T,Mtyp,rhom,cpm0);  % organic-rich material
 rhocpi = C_ice(T);                     % ice
 rhocpu = C_water(T);                   % unfrozen water

% find the lattice-vibration component of the heat capacity

 rhocp = (1 - phi).*rhocpm + (phi_i).*rhocpi + (phi_u).*rhocpu;

% find the phase change contribution

 Cphase = rhow * Lf * dphiudT;

% find the total volumetric heat capacity

 C = rhocp  + Cphase;
