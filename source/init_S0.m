 function S0 = init_S0(Marr,col,CS)

% Initializes the heat-source parameter S0.
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

%   Marr = array of material properties                         (M+1,nparams)
%   col  = column within Marr where So is located               (scalar)
%   CS   = coordinate system                                    (character string)
%   S0   = heat-source extrapolated to the surface              (M+1)
% ______________________________________________

 switch CS
 case 'R'
   Marr = Marr';
 end

% extract S0 from the Marr array

 S0 = Marr(:,col);

% convert to SI base units

 S0 = 1e-06 * S0;       % [W/m^3]
