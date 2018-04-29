 function cpm0 = init_cpm(Mtyp,Marr,col,CS)

% Initializes the specific heat of the matrix material at reference 
% temperature (20 C).
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
%   col  = column within Marr where cpm0 is located (scalar)
%   CS   = coordinate system                        (character string)
%   cpm0 = matrix specific heat at 20 C             (M+1)

% Notes:
%   (1) For ice (Mtyp=3), ignore what's found in the *_layers.txt file.
% ______________________________________________

% set cpm0 according to what's found in the *_layers.txt file

 switch CS
 case 'R'
   cpm0 = Marr(col,:);
 otherwise
   cpm0 = Marr(:,col);
 end

% reset cmp0 if layer consists of pure ice

 L = Mtyp == 3;

 if any(L)
   cpm0(L) = NaN;
 end
