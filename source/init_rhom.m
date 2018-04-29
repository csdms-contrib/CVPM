 function rhom = init_rhom(Mtyp,Marr,col,CS)

% Initializes the density of the matrix material.
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
%   Marr = array of material parameters             (M+1,nparams)
%   col  = column within Marr where rhm0 is located (scalar)
%   CS   = coordinate system                        (character string)
%   rhom = matrix density                           (M+1)

% Notes:
%   (1) For ice (Mtyp=3), ignore what's found in the *_layers.txt file.
% ______________________________________________

% set rhom according to what's found in the *_layers.txt file

 switch CS
 case 'R'
   rhom = Marr(col,:);
 otherwise
   rhom = Marr(:,col);
 end

% reset rhom if layer consists of pure ice

 L = (Mtyp >= 2) & (Mtyp <= 3);

 if any(L)
   rhom(L) = 917;
 end
