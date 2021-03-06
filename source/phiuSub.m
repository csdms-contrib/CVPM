 function phi_u = phiuSub(Mtyp,phi_w,Psi1,Psi2,r1,r2,lambda,DeltaT)

% Finds the volume fraction of unfrozen water during the initial setup
% by CPS.
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

%   Mtyp    = material type                                     (M+1)
%   phi_w   = volume fraction of water (solid + liquid)         (M+1)
%   Psi1    = relative volume fraction, large particles & pores (M+1)
%   Psi2    = relative volume fraction, small particles & pores (M+1)
%   r1      = effective radius of large particles or pores      (M+1)
%   r2      = effective radius of small particles or pores      (M+1)
%   lambda  = interfacial melting parameter                     (M+1)
%   DeltaT  = Tf - T                                            (M+1)
%   phi_u   = volume fraction of unfrozen water                 (M+1)
% ______________________________________________

% parameters

 PC_opt = 'fcc';

% pre-allocate arrays

 phi_u  = zeros(size(lambda));
 phi_u1 = zeros(size(lambda));
 phi_u2 = zeros(size(lambda));

% find which material types are present

 [Mtypes,m] = which_Mtypes(Mtyp);

% find the volume fraction of water for each particle or pore size

 phi_w1 = Psi1 .* phi_w;
 phi_w2 = Psi2 .* phi_w;

% > Find volume fraction of unfrozen water

 for j=1:m      % loop through unique material types

   L = Mtyp == Mtypes(j);

   switch fix(Mtypes(j))
   case 1                           % fixed properties (used for testing)
     phi_u(L) = 0;

   case {2,3}                       % pure ice
     phi_u(L) = 0;

   case {4,5,6,7,8,9,10,11,12,13,14,15,20,21,22,23,24,25}   % rocks and organic-rich materials

     phi_u1(L) = phiu_pore(phi_w1(L),r1(L),lambda(L),DeltaT(L),PC_opt);
     phi_u2(L) = phiu_pore(phi_w2(L),r2(L),lambda(L),DeltaT(L),PC_opt);
     phi_u(L)  = phi_u1(L) + phi_u2(L);
   end
 end
