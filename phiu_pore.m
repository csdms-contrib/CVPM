 function phi_u = phiu_pore(phi_w,r,lambda,DeltaT,PC_opt)

% Finds the volume fraction of unfrozen water in pores for matrix particles
% with an effective radius r based on the Cahn et al. (1992) model.
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

%   phi_w  = volume fraction of water (solid & liquid)  (vector)
%   r      = effective radius of matrix particles [m]   (vector)
%   lambda = interfacial melting parameter [K^1/3 m]    (vector)
%   DeltaT = Tf - T                                     (vector)
%   PC_opt = packing coefficient option                 (string)
%   phi_u  = volume fraction of unfrozen water          (vector)

% Notes:

%   PC_opt = 'sc' or 'fcc'
% ______________________________________________

% parameters

 xi = 0.0259e-06;   % curvature coefficient [K m]

 switch PC_opt
 case 'sc'          % simple cubic packing (sc)
   a1 = 1.893;
   a2 = 3.367;
 case 'fcc'         % cubic close packing (fcc)
   a1 = 2.450;
   a2 = 8.572;
 otherwise
   disp(' ')
   disp('warning phiu_pore: unknown packing option.')
   pause
 end

% > Find phi_u

% reset DeltaT values corresponding to temperatures warmer than Tf to yield
% phi_u values equal to phi 

 L         = DeltaT <= 0;
 DeltaT(L) = 1e-09;

% term1: water film on surface of spheres and at grain boundaries

 t1    = lambda ./ (r .* DeltaT.^(1/3));
 phiu1 = a1 * t1;

% term2: water in crevices and grain-boundary edges

 t2    = xi ./ (r .* DeltaT);
 phiu2 = a2 * t2.^2;

% add terms

 phi_u = phiu1 + phiu2;

% phi_u cannot exceed phi_w

 L = phi_u > phi_w;        
 if length(phi_w) == 1
   phi_u(L) = phi_w;
 elseif length(phi_w) == length(phi_u)
   phi_u(L) = phi_w(L);
 else
   disp(' ')
   disp('phiu_pore: warning, unknown case.')
   pause
 end

% reset if r = 0

% L = r == 0;
% phi_u(L) = 0;
