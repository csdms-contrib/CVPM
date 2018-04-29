 function aw = water_activity(solute,xs)

% Finds the activity of water (aw) for aqueous solutions with various
% electrolytes using the expression provided by Miyawaki et al (1997). 
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

% To prevent the code from crashing, aw values are linearly extrapolated
% when solute concentrations are greater than the max(xs) value for which
% the Miyawaki coefficients are valid.

% Notation:

%   solute  = chemical formula for the solute   (string)
%   xs      = solute mole fraction              (vector)
%   aw      = water activity                    (vector)

% Currently available solutes include,

%   'NaCl'
%   'KCl'
% ______________________________________________

 switch solute
 case 'NaCl'
   alpha  =  1.825;
   beta   = -20.78;
   xsmax  =  0.15;      % max(xs) value for alpha,beta validity
   aw_xs  = 0.82565;    % aw @ max(xs)
   dawdxs = -1.6775;    % daw/dxs at max(xs)
 case 'KCl'
   alpha  =  4.754;
   beta   = -49.37;
   xsmax  =  0.07;
   aw_xs  = 0.93593;
   dawdxs = -1.0637;
 otherwise
   disp(' ')
   disp('Water_activity: this solute is currently unavailable.')
   disp(' ')
   pause
 end

% Miyawaki equation

 aw = (1 - xs) .* exp(alpha * xs.^2 + beta * xs.^3);

% use linear extrapolation beyond xsmax

 L     = xs > xsmax;
 aw(L) = aw_xs + dawdxs*(xs(L) - xsmax);

% aw can't be negative

 L     = aw < 0;
 aw(L) = 0;
