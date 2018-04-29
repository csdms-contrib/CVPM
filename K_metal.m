 function K = K_metal(T,Mtyp)

% Finds the thermal conductivity of metals.
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

%   T    = temperature (C)          (vector)
%   Mtyp = metal type               (scalar)
%   K    = thermal conductivity     (vector)

% Available metals:

%   40  steel drill pipe
%   41  stainless steel
%   42  cast iron
%   43  aluminum
%   44  copper

% Notes:
%   (1) This routine does not currently take into account the temperature 
%       dependence.

% Source: wikipedia (for all values except steel drill pipe)
% ______________________________________________

 switch Mtyp(1)
 case 40            % steel drill pipe
   K = 45;

 case 41            % stainless steel (note: values vary considerably)
   K = 16.5;

 case 42            % cast iron
   K = 55;

 case 43            % aluminum
   K = 236;

 case 44            % copper
   K = 401;
 end

 K = K * ones(size(Mtyp));
