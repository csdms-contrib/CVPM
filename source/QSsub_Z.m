 function QS = QSsub_Z(S_opt,S0,hs,zf)

% Initializes the integral of the heat source over the control volumes
% (vertical dimension).
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

% Currently available heat sources:

%   S_opt = 'zero'           no heating
%         = 'linear'         source decreases linearly with depth
%         = 'exponential'    source decreases exponentially with depth

% To allow for a variable lithologic profile, different layers are allowed to 
% have different S0 and hs values.

% Notation:

%   S_opt = source function option                  (character string)
%   S0    = source function at Z = 0                (M+1)
%   hs    = source function length scale            (M+1)
%   zf    = CV interfaces                           (M+1)
%   QS    = source function integral                (M+1)
% ______________________________________________

 M = length(zf) - 1;

% pre-allocate array

 QS = NaN*ones(size(zf));

% integrate S(z) over the control volumes

 switch S_opt
 case 'zero'
   QS = zeros(size(zf));

 case 'linear'
   for k=2:M
     QS(k) = S0(k) * ((zf(k+1) - zf(k)) - (zf(k+1)^2 - zf(k)^2)/(2*hs(k)));
   end

 case 'exponential'
   for k=2:M
     QS(k) = S0(k)*hs(k) * (expm1(-zf(k)/hs(k)) - expm1(-zf(k+1)/hs(k)));
   end
 end
