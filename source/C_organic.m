 function rhocp = C_organic(T,Mtyp,rhom,cpm0)

% Finds the volumetric heat capacity of a organic-rich matrix at temperature
% T.  The matrix is assumed to consist of a mixture of peat and mineral
% particles.
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

%   T     = temperature (C)                     (vector)
%   Mtyp  = material type                       (vector)
%   rhom  = mean density of mineral grains      (vector)
%   cpm0  = specific heat of minerals (20 C)    (vector)
%   rhocp = volumetric heat capacity            (vector)

% Available matrix types (Mtyp):

%   20      100% peat
%   21  	 75% peat/25% mineral mix
%   22       50% peat/50% mineral mix
%   23       25% peat/75% mineral mix
% ______________________________________________

 Tk = T + 273.15;       % kelvin temperature

% heat capacity of peat particles (see 2004Ling_a).

 rhocpo = 1e06*(-0.1333 + 6.255e-03 * Tk);

% heat capacity of mineral particles

 rhocpm = C_mineral(T,Mtyp,rhom,cpm0);

% heat capacity of mixture

 switch Mtyp(1)
 case 20
   f = 1;
 case 21
   f = 0.75;
 case 22
   f = 0.5;
 case 23
   f = 0.25;
 case 24
   f = 0;
 end

 rhocp = f *rhocpo + (1-f) *rhocpm;
