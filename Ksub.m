 function K = Ksub(T,Mtyp,Km0,phi,phi_i,phi_u,planet)

% Finds the bulk thermal conductivity for each CV.  For convenience, K
% is also defined on the boundaries where it is set equal to that of the 
% adjacent CVs.
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

%   T     = temperature (C)                         (M+1)
%   Mtyp  = material type                           (M+1)
%   Km0   = matrix conductivity at 0 C              (M+1)
%   phi   = porosity                                (M+1)
%   phi_i = volume fraction of ice                  (M+1)
%   phi_u = volume fraction of unfrozen water       (M+1)
%   K     = bulk thermal conductivity               (M+1)
% ______________________________________________

 M = length(T) - 1;

% pre-allocate arrays

 K = NaN*ones(size(T));

% find which material types are present

 [Mtypes,m] = which_Mtypes(Mtyp);

% find K for each material present

 for j=1:m
   L = Mtyp == Mtypes(j);

   switch fix(Mtypes(j))
   case 1                       % fixed properties (used for testing)

     K(L) = Km0(L);

   case 2                       % linerized ice (used for testing)

     K(L) = K_linIce(T(L),Km0(L));

   case 3                       % ice

     K(L) = K_ice(T(L));

   case {4,5,6,7,8,9,10,11,12,13,14,15} % rocks and soils

     K(L) = K_rock(T(L),Mtyp(L),Km0(L),phi(L),phi_i(L),phi_u(L),planet);

   case {20,21,22,23,24,25}     % organic-rich material

     K(L) = K_veg(T(L),Mtyp(L),Km0(L),phi(L),phi_i(L),phi_u(L),planet);

   case {30,31,32,33,34,35}     % borehole fluids

     Ftype = Mtypes(j) - 29;
     K(L)  = K_fluid(T(L),Ftype);

   case {40,41,42,43,44}        % metals

     K(L) = K_metal(T(L),Mtyp(L));
   end
 end
