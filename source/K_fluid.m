 function K = K_fluid(T,Ftype)

% Finds the thermal conductivity of borehole fluids.
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

%   T     = temperature (C)          (vector)
%   Ftype = fluid type               (scalar)
%   K     = thermal conductivity     (vector)

% Available fluids:

%   1  water
%   2  JetA
%   3  n-butyl acetate
%   4  Estisol 140
%   5  Estisol 240
%   6  IsoparK
% ______________________________________________

 switch Ftype
 case 1                     % water
   K = K_water(T);

 case 2                     % JetA (at -10 C)
   K = 0.12046*ones(size(T));

 case 3                     % n-butyl acetate (at -10 C)
   K  = 0.1434*ones(size(T));

 case {4,5}                 % Estisol 140/240
   K = 0.18*ones(size(T));

 case 6                     % IsoparK mix
   K = 0.15*ones(size(T));
 end
