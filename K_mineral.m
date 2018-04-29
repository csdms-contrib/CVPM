 function Km = K_mineral(T,Mtyp,Km0)

% Finds the thermal conductivity of mineral grains at temperature T.
% This routine takes into account the difference between 
% igneous/sedimentary and sedimentary mineral assemblages.
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

%   T    = temperature (C)                              (vector)
%   Mtyp = material type                                (vector)
%   Km0  = mineral conductivity at 0 C                  (vector)
%   Km   = mineral conductivity at temperature T        (vector)
% ______________________________________________

% pre-allocate the array

 Km = NaN*ones(size(Km0));

% igneous/metamorphic minerals

 L = (Mtyp >= 4) & (Mtyp <= 9);

 if any(L)
   a = [0.99 0.0030 0.0042];
   Km(L) = Km0(L) ./ (a(1) + T(L) .* (a(2) - a(3)./Km0(L)));
 end

% sedimentary minerals

 L = (Mtyp >= 10) & (Mtyp <= 15);

 if any(L)
   a = [0.99 0.0034 0.0039];
   Km(L) = Km0(L) ./ (a(1) + T(L) .* (a(2) - a(3)./Km0(L)));
 end
