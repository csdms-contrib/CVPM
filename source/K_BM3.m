 function K = K_BM3(K1,K2,K3,psi1,psi2,psi3)

% Finds the bulk thermal conductivity of a 3-component system using the 
% Brailsford and Major model that assumes phases 2 and 3 are randomly
% distributed within a continuous phase 1 medium.
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

%   K1   = conductivity of phase 1                 (vector)
%   K2   = conductivity of phase 2                 (vector)
%   K3   = conductivity of phase 3                 (vector)
%   psi1 = relative volume fraction of phase 1     (vector)
%   psi2 = relative volume fraction of phase 2     (vector)
%   psi3 = relative volume fraction of phase 3     (vector)
%   K    = bulk thermal conductivity               (vector)

% Notes:
%   (1)   1 = psi1 + psi2 + psi3      (this should be true)
% ______________________________________________

 chi2 = (K1 ./ K2);
 chi3 = (K1 ./ K3);

 f2   = psi2 ./ (2*chi2 + 1);
 f3   = psi3 ./ (2*chi3 + 1);
 f4   = f2 + f3;
 f5   = f2 .* chi2 + f3 .* chi3;

 K    = K1 .* (psi1 + 3 * f4) ./ (psi1 + 3 * f5);
