 function r = init_r(Mtyp,Marr,col,CS)

% Initializes the effective radius of the matrix particles.
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
%   col  = column within Marr where r is located	(scalar)
%   CS   = coordinate system                        (character string)
%   r    = particle radius                          (M+1)
% ______________________________________________
 
 switch CS
 case 'R'
   Marr = Marr';
 end

% pre-allocate diameters assuming no pores are present

 d = zeros(size(Mtyp));

% reset the radius for rocks, soils, and organic-rich materials

 L = (Mtyp >= 4) & (Mtyp <= 25);

 if any(L)
   d(L) = Marr(L,col);
 end

 r = d/2;

% convert to SI base units

 r = 1e-06 * r;     % [m]
