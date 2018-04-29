 function theta_p = dT_press(planet,P_opt,Z,dz,rho)

% Finds the freezing point depression due to pressure.  
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

% Pore pressures can be calculated assuming either 'hydrostatic' or 
% 'lithostatic' environmental conditions.

% Notation:

%   planet                                   (string)
%   P_opt   = pressure option                (string)
%   Z       = depth                          (M+1)
%   dz      = thickness of each CV           (M+1)
%   rho     = bulk density                   (M+1)
%   theta_p = feezing point depression       (M+1)

% Notes:

% (1) planet = 'earth' or 'mars'
% (2) P_opt  = 'hydrostatic' or 'lithostatic'
% ______________________________________________

% parameters

 a    = 9.8e-08;        % air-saturated pressure sensitivity (K/Pa)
 P0   = 611.657;        % triple point pressure (Pa)
 rhow = 1000;           % density of water (kg/m^3)

 switch planet
 case 'earth'
   g    = 9.78;         % surface gravity (m/s^2)
   Patm = 1.01e05;      % surface air pressure (Pa)
 case 'mars'
   g    = 3.72;
   Patm = 700;
 otherwise
   disp(' ')
   disp('warning dT_press: unknown planet.')
   pause
 end

% find pore pressure

 switch P_opt
 case 'off'
   P = P0 *ones(size(Z));
 case 'hydrostatic'
   rhoZ = rhow * Z;
   P    = Patm + g * rhoZ;
 case 'lithostatic'
   rhodz = rho .* dz;
   rhoZ  = zeros(size(Z));
   for i=1:length(Z)
     rhoZ(i) = sum(rhodz(1:i));
   end
   P = Patm + g * rhoZ;
 otherwise
   disp(' ')
   disp('warning dT_press: unknown pressure option.')
   pause
 end

% find freezing point depression

 theta_p = a * (P - P0);
