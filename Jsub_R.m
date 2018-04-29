 function [Jtilde,Jnet] = Jsub_R(BCtypes,T,rf,Ke,dr,qo)

% Finds the integrated heat fluxes across CV interfaces and the net
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

% J and Jnet calculated by this function are not used to find the temperture
% field.  Instead, they are provided to assist with CVPM diagnostics.

% Notation:

%   T      = temperature                              (N+1)
%   rf     = interface locations                      (N+1)
%   Ke     = effective conduct. across R-interfaces   (N+1)
%   qo     = heat flux through outer interface        (scalar)
%   Jtilde = integrated heat flux at CV interfaces    (N+1)
%   Jnet   = net heat-flux out of each CV             (N+1)

% Notes: 
%   (1) qb is not used if lowerBC_type = 'T'.
%   (2) For the inner CV, the heat flux J(2) initially calculated 
%       from the temperature difference T(2)-T(1) pertains to a point
%       midway between grid points W and P.  This value and J(3) are
%       are then used to find the flux at the inner boundary.  This
%       is then stored in J(2).
% ______________________________________________

 N = length(T) - 1;
 innerBC_type = BCtypes{1};
 outerBC_type = BCtypes{2};

% pre-allocate arrays

 Jtilde = NaN*ones(size(T));
 Jnet   = NaN*ones(size(T));

% find the heat flux across each CV interface

 rKedr = rf .* Ke ./ dr;

 for k=2:N
   Jtilde(k) = -rKedr(k) * (T(k) - T(k-1));
 end

% find an improved value for the flux on the inner boundary

 switch innerBC_type
 case 'T'
   Jtilde(2) = (4*Jtilde(2) - Jtilde(3)) / 3;
 end

 switch outerBC_type
 case 'q'
   Jtilde(N+1) = 2*pi*rf(N+1)*qo;
 end

% find net heat-flux out of each CV

 for k=2:N
   Jnet(k) = Jtilde(k+1) - Jtilde(k);
 end
