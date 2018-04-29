 function [KeR,KeZ] = Keff_RZ(K,varepR,varepZ)

% Finds the effective conductivity at the CV interfaces for the RZ
% coordinate system.
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

%   K      = conductivity at CV grid points                 (M+1,N+1)
%   varepR = fractional distance radial direction           (N+1)
%   varepZ = fractional distance vertical direction         (M+1)
%   KeR    = effective conductivity across R-interfaces     (M+1,N+1)
%   KeZ    = effective conductivity across Z-interfaces     (M+1,N+1)
% ______________________________________________

 N = length(varepR) - 1;
 M = length(varepZ) - 1;
 
% pre-allocate array

 KeR = NaN*ones(size(K));
 KeZ = NaN*ones(size(K));

% find effective conductivity at R-interfaces

 for k=1:M                      % do one row (k) at a time
   KR       = K(k,:);
   KeR(k,:) = Keff(KR,varepR);
 end
 KeR(M+1,:) = KeR(M,:);

% find effective conductivity at Z-interfaces

 for j=1:N                      % do one column (j) at a time
   KZ       = K(:,j);
   KeZ(:,j) = Keff(KZ,varepZ);
 end
 KeZ(:,N+1) = KeZ(:,N);
