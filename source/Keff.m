 function Ke = Keff(K,varep)

% Finds the effective conductivity at CV interfaces.
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

%   K     = conductivity at the CV grid pts             (M+1)
%   varep = fractional distance                         (M+1)
%   Ke    = effective conductivity at CV interfaces     (M+1)
% ______________________________________________

 M = length(K) - 1;

% pre-allocate array

 Ke = NaN*ones(size(varep));

% find effective conductivity

 for k=3:M
   Ke(k) = 1 / ( (1-varep(k))/K(k-1) + varep(k)/K(k) );
 end

 Ke(1:2) = K(2);    % effective Ke on upper boundary
 Ke(M+1) = K(M);    %     "     "     lower boundary
