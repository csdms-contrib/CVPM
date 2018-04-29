 function K = K_veg(T,Mtyp,Km0,phi,phi_i,phi_u,planet)

% Finds the bulk thermal conductivity of organic-rich material at 
% temperature T.  The material is assumed to consist of a matrix of 
% organic and mineral particles while the pores are filled with liquid 
% water, ice, and air.
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
%   Mtyp  = material type                           (vector)
%   Km0   = mineral conductivity at 0 C             (vector)
%   phi   = porosity                                (vector)
%   phi_i = volume fraction of ice                  (vector)
%   phi_u = volume fraction of unfrozen water       (vector)
%   K     = bulk thermal conductivity               (vector)

% Note:
%   (1) This routine is identical to K_rock except the matrix has both
%       organic and mineral particles.
% ______________________________________________

% pre-allocate arrays

 K = NaN*ones(size(phi));

% find the conductivity of the organic matrix

 Km = K_organic(T,Mtyp,Km0);

% find the bulk conductivity of porous material (use same code as for soils)

 K = K_rock_porous(T,Km,phi,phi_i,phi_u,planet);
