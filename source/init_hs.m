 function hs = init_hs(Marr,col,CS)

% Initializes the heat-source length scale.
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

%   Marr = array of material properties             (M+1,nparams)
%   col  = column within Marr where So is located	(scalar)
%   CS   = coordinate system                        (character string)
%   hs   = heat-source length scale              	(M+1)
% ______________________________________________

 switch CS
 case 'R'
   Marr = Marr';
 end

% extract hs from the Marr array

 hs = Marr(:,col);

% convert to SI base units

 hs = 1000 * hs;        % [m]
