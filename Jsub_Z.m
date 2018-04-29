 function [J,Jnet] = Jsub_Z(BCtypes,T,Ke,dz,qb)

% Finds the diffusive heat fluxes across CV interfaces and the net
% heat flux out of each control volume.  
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

% J and Jnet calculated by this function are not used to find the temperature
% field.  Instead, they are provided to assist with CVPM diagnostics.

% Notation:

%   T    = temperature                              (M+1)
%   Ke   = effective conductivity                   (M+1)
%   dz   = distance between grid points             (M+1)
%   qb   = heat flux through bottom interface       (scalar)
%   J    = diffusive heat flux at CV interfaces     (M+1)
%   Jnet = net heat-flux out of each CV             (M+1)

% Notes: 
%   (1) qb is not used if lowerBC_type = 'T'.
%   (2) For the upper CV, the heat flux J(2) initially calculated 
%       from the temperature difference T(2)-T(1) pertains to a point
%       midway between grid points U and P.  This value and J(3) are
%       are then used to find the flux at the upper boundary.  This
%       is then stored in J(2).
% ______________________________________________

 M            = length(T) - 1;
 Kedz         = Ke ./ dz;
 upperBC_type = BCtypes{1};
 lowerBC_type = BCtypes{2};

% pre-allocate arrays

 J    = NaN*ones(size(T));
 Jnet = NaN*ones(size(T));

% find the heat flux across each CV interface

 for k=2:M
   J(k) = -Kedz(k) * (T(k) - T(k-1));
 end

% find an improved value for the flux on the upper boundary

 switch upperBC_type
 case 'T'
   J(2) = (4*J(2) - J(3)) / 3;
 end

 switch lowerBC_type
 case 'q'
   J(M+1) = qb;
 end

% find net heat-flux out of each CV

 for k=2:M
   Jnet(k) = J(k+1) - J(k);
 end
