 function Km0 = init_Km(Mtyp,Marr,col,CS)

% Initializes the thermal conductivity of the matrix material at reference 
% temperature 0 C.
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
%   col  = column within Marr where Km0 is located  (scalar)
%   CS   = coordinate system                        (character string)
%   Km0  = matrix conductivity at 0 C               (M+1)

% Notes:
%   (1) The contents of Marr are ignored when setting Km0 for layers
%       consisting of ice or borehole fluid.
% ______________________________________________

% set Km0 according to what's found in the *_layers.* file

 switch CS
 case 'R'
   Km0 = Marr(col,:);
 otherwise
   Km0 = Marr(:,col);
 end

% reset Km0 if layer consists of pure ice

 L = Mtyp == 3;

 if any(L)
   Km0(L) = K_ice(0);
 end

% reset Km0 if layer consists of borehole fluid

 L = (Mtyp >= 30) & (Mtyp <= 35);

 if any(L)
   Ftyp   = Mtyp(L) - 29;
   Ftype  = unique(Ftyp);
   Km0(L) = K_fluid(0,Ftype);
 end
