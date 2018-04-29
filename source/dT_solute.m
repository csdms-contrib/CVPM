 function theta_s = dT_solute(Mtyp,solute,xs0,phi_w,phi_u)

% Finds the freezing point depression due to solutes.
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

%   Mtyp    = material type                             (M+1)
%   solute  = chemical formula for the solute           (string)
%   xs0     = solute mole fraction with no ice          (M+1)
%   phi_w   = volume fraction of water (solid+liquid)   (M+1)
%   phi_u   = volume fraction of unfrozen water         (M+1)
%   theta_s = freezing point depression                 (M+1)

 dflag = 0;     % (0) quiet mode, (1) dumps various fields to screen.
% ______________________________________________

% pre-allocate arrays

 theta_s = zeros(size(Mtyp));

% find the solute mole fraction

 xs = xs0 .* (phi_w ./ phi_u);

% find the freezing point depression

 L = (Mtyp >=4) & (Mtyp <= 25);

 if any(L)
   theta_s(L) = dT_solute_sub(solute,xs(L));
 end

 if dflag
   disp(' ')
   disp('dT_solute: Mtyp phi_w phi_u xs theta_s')
   [Mtyp phi_w phi_u xs theta_s]
   pause(2)
 end
