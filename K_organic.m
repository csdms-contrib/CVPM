 function K = K_organic(T,Mtyp,Km0)

% Finds the thermal conductivity of an organic-rich matrix at temperature T.
% The matrix is assumed to consist of a mixture of peat and mineral particles.
% The mineral particles are assumed to be those typically found in mud rocks
% (mineral Mtyp=11).
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

% Available matrix types:

%   20      100% peat
%   21  	 75% peat/25% mineral mix
%   22       50% peat/50% mineral mix
%   23       25% peat/75% mineral mix
% ______________________________________________

% find the conductivity of the peat particles (see 2004Ling_a)

 Ko = 0.25 *ones(size(Km0));

% find the conductivity of the mineral particles

 mtyp = 11 *ones(size(Km0));
 Km   = K_mineral(T,mtyp,Km0);

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
 f = f *ones(size(Km0));

% conductivity of a random mixture

 K = K_BM2r(Km,Ko,f);
