 function [phi_u,T] = phiu_table(Mtyp,phi_w,Psi1,Psi2,r1,r2,lambda,solute,xs0,theta_p)

% Builds a table of phi_u values corresponding to a DeltaT array (ranging 
% from 500 C to -100 C) at each CV.  
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

% The temperatures T corresponding to the DeltaT array are then determined 
% (this will be different for each layer due to the pressure and solute 
% effects).

% CPS subsequently interpolates the values in this table to find
% the unfrozen water content (phi_u) at temperature T.

% Notation:

%   Mtyp    = material type                                     (M+1)
%   phi_w   = volume fraction of water (solid & liquid)         (M+1)
%   Psi1    = relative volume fraction, large particles & pores (M+1)
%   Psi2    = relative volume fraction, small particles & pores (M+1)
%   r1      = effective radius of large particles or pores      (M+1)
%   r2      = effective radius of small particles or pores      (M+1)
%   lambda  = interfacial melting parameter                     (M+1)
%   solute  = chemical formula for the solute                   (string)
%   xs0     = solute mole fraction with no ice                  (M+1)
%   theta_p = pressure freezing-point depression                (M+1)
%   phi_u   = volume fraction of unfrozen water                 (M+1,nT)
%   T       = temperatures corresponding to phi_u               (M+1,nT)

% Notes:

%   (1) The volume fraction of water (phi_w) is constrained by the porosity
%       (phi) and the degree of saturation (Sr) at the time of the initial
%       condition.  As time proceeds, phi_w will change in response to 
%       changing temperatures as ice converts to unfrozen water and vice 
%       versa.  As phi_w decreases upon warming and increases upon cooling,
%       phi_u can never exceed the initial phi_w value.  Thus, this table is
%       constructed within CPS using the constraint phi_u <= initial(phi_w).
% ______________________________________________

 Mp1 = length(Mtyp);    % M+1

% parameters

 Tf0    = 0.01;     % triple point temperature (C)
 PC_opt = 'fcc';    % packing coefficient option

% define DeltaT array

 x      = [-500 -100 -10 -1 -0.1 0];        % 500 to 0 C
 k      = (-3:0.005:2);
 y      = 10.^k;                            % -0.001 to -100 C
 DeltaT = [x y];
 nT     = length(DeltaT);

% find the volume fraction of water for each particle or pore size

 phi_w1 = Psi1 .* phi_w;
 phi_w2 = Psi2 .* phi_w;

% pre-allocate arrays

 phi_u  = zeros(Mp1,nT);
 phi_u1 = zeros(Mp1,nT);
 phi_u2 = zeros(Mp1,nT);
 T      = NaN*ones(size(phi_u));

% Step through layers

 for k=1:Mp1       % step through layers

   if (Mtyp(k) >= 4) && (Mtyp(k) <= 25)

%   find phi_u values corresponding to DeltaT array

     phi_u1(k,:) = phiu_pore(phi_w1(k),r1(k),lambda(k),DeltaT,PC_opt);
     phi_u2(k,:) = phiu_pore(phi_w2(k),r2(k),lambda(k),DeltaT,PC_opt);
     phi_u( k,:) = phi_u1(k,:) + phi_u2(k,:);

%   find temperatures T corresponding to the DeltaT array

     xs      = xs0(k) .* (phi_w(k) ./ phi_u(k,:));  % solute mole fraction at DeltaTA
     theta_s = dT_solute_sub(solute,xs);            % solute freezing pt depression at DeltaTA
     Tf      = Tf0 - (theta_p(k) + theta_s);        % bulk freezing temperature at DeltaTA
     T(k,:)  = Tf - DeltaT;                         % temperatures corresponding to DeltaTA
   end
 end
