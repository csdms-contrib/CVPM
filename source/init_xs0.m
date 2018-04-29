 function xs0 = init_xs0(Mtyp,Marr,col,CS)

% Initializes the mole fraction of any solutes in the pore spaces if all
% of the pore ice were to melt.
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

%   Mtyp = material type                            (M+1)
%   Marr = array of material properties             (M+1,nparams)
%   col  = column within Marr where xs0 is located  (scalar)
%   CS   = coordinate system                        (character string)
%   xs0  = mole fraction of solutes with no ice     (M+1)
% ______________________________________________

 switch CS
 case 'R'
   Marr = Marr';
 end

% pre-allocate the array assuming no solutes are present

 xs0 = zeros(size(Mtyp));

% reset the mole fraction for rocks, soils, and organic-rich materials

 L = (Mtyp >= 4) & (Mtyp <= 25);

 if any(L)
   xs0(L) = Marr(L,col);
 end
