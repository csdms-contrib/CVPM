 function theta_s = dT_solute_sub(solute,xs)

% Finds the freezing point depression due to solutes for non-ideal
% solutions using the expression reported by Robinson and Stokes (1959).
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

%   solute  = chemical formula for the solute   (string)
%   xs      = solute mole fraction              (M+1)
%   theta_s = freezing point depression         (M+1)
% ______________________________________________

% parameters

 a = 4.76e-06;
 b = 9.687e-3;

% if solutes are present ...

 if ~strcmp(solute,'none')

   switch solute
   case 'NaCl'
     T_eutectic = -21.1;
   otherwise
     disp(' ')
     disp('warning dT_solute: eutectic temperature for this solute is unknown.')
     pause
   end

% find the water activity

   aw = water_activity(solute,xs);
   c  = log(aw);

% find the freezing point depression

   theta_s = (-b + sqrt(b^2 - 4*a*c))/(2*a);

% check for eutectic pt

   L          = theta_s > -T_eutectic;
   theta_s(L) = -T_eutectic;

 else
   theta_s = zeros(size(xs));
 end
