 function K = K_BM2r(K1,K2,psi2)

% Finds the bulk thermal conductivity of a 2-component system using the 
% Brailsford and Major model that assumes both phases are randomly
% distributed (their eq. 13, corrected).
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

%   K1   = conductivity of phase 1                  (vector)
%   K2   = conductivity of phase 2                  (vector)
%   psi2 = relative volume fraction of phase 2      (vector)
%   K    = bulk thermal conductivity                (vector)

% Notes:
%   (1)     1 = psi1 + psi2
% ______________________________________________

 chi = (K1 ./ K2);

 f1  = 2*chi - 1;
 f2  = 3*psi2 .* (chi - 1);
 f   = f1 - f2;

 K   = K1 .* (f + sqrt(f.^2 + 8*chi)) ./ (4*chi);
