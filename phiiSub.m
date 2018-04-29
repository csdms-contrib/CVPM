 function phi_i = phiiSub(Mtyp,Mw,phi_u)

% Finds the volume fraction of ice once Mw is known.
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

%   Mtyp    = material type                                 (M+1)
%   Mw      = mass of water (both phases) per unit volume   (M+1)
%   phi_u   = volume fraction of unfrozen water             (M+1)
%   phi_i   = volume fraction of ice                        (M+1)
% ______________________________________________

% parameters

 rhoi = 917;            % density of ice
 rhou = 1000;           % density of water

% pre-allocate arrays

 phi_i = zeros(size(phi_u));

% find which material types are present

 [Mtypes,m] = which_Mtypes(Mtyp);

% > Find volume fraction of ice

 for j=1:m      % loop through unique material types

   L = Mtyp == Mtypes(j);

   switch fix(Mtypes(j))
   case 1                           % fixed properties (used for testing)
     phi_i(L) = 0;

   case {2,3}                       % pure ice
     phi_i(L) = 1;

   case {4,5,6,7,8,9,10,11,12,13,14,15,20,21,22,23,24,25}   % rocks and organic-rich materials

     phi_i(L)  = (Mw(L) - rhou *phi_u(L)) / rhoi;
     LL        = phi_i < 0;
     phi_i(LL) = 0;
   end
 end
