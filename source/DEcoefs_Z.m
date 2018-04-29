 function [aU,aUp,aD,aDp,aP,aPp,b] = DEcoefs_Z(dt,f,Az,VP,C,Ke,QS,BCtypes,q,qn)

% Finds the discretization coefficients for 1-D vertical CVPM model.
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

% Upper boundary condition is set by,

%   upperBC_type:   'T' prescribed temperature
%                   'q' prescribed heat flux

% Lower boundary condition is set by,

%   lowerBC_type:   'T' prescribed temperature
%                   'q' prescribed heat flux

% Notation:

%   dt                      time step                           (scalar)
%   f                       implicit/explicit factor            (scalar)
%   Az                      1/dz                                (M+1)
%   VP                      Dz                                  (M+1)
%   C                       volumetric heat capacity            (M+1)
%   Ke                      effective conduct., Z-interfaces    (M+1)
%   QS                      spatially integrated source term    (M+1)
%   BCtypes                 upper/lower BC types                (string)
%   q                       heat fluxes, this time step         (2)
%   qn                      heat fluxes, last time step         (2)
%   aU,aUp,aD,aDp,aP,aPp,b  DE coefficients                     (M+1)
% ______________________________________________

 M = length(C) - 1;

% pre-allocate arrays

 aU  = zeros(M+1,1);
 aUp = zeros(M+1,1);
 aD  = zeros(M+1,1);
 aDp = zeros(M+1,1);

% unpack BCs

 upperBC_type = BCtypes{1};
 lowerBC_type = BCtypes{2};

 switch upperBC_type
 case 'q'
   qs  = q(1);
   qsn = qn(1);
 end
 switch lowerBC_type
 case 'q'
   qb  = q(2);
   qbn = qn(2);
 end

% indices of upper/lower CVs

 kU = 2;        % upper CV
 kL = M;        % lower CV

% define the following for convenience

 g     = 4/3;
 dtf   = dt * f;
 dt1mf = dt * (1-f);

 VPC   = VP .* C;
 AKeZ  = Az .* Ke;

% coefficients aU,aUp,aD,aDp

 aU( kU:kL) = dtf   * AKeZ(kU:kL);
 aUp(kU:kL) = dt1mf * AKeZ(kU:kL);
 aD( kU:kL) = dtf   * AKeZ(kU+1:kL+1);
 aDp(kU:kL) = dt1mf * AKeZ(kU+1:kL+1);

 k = kU;
 switch upperBC_type
 case 'T'
   aU( k) = g * aU( k);
   aUp(k) = g * aUp(k);
   aD( k) = g * aD( k);
   aDp(k) = g * aDp(k);
 case 'q'
   aU( k) = 0;
   aUp(k) = 0;
 end

 k = kL;
 switch lowerBC_type
 case 'T'
   aU( k) = g * aU( k);
   aUp(k) = g * aUp(k);
   aD( k) = g * aD( k);
   aDp(k) = g * aDp(k);
 case 'q'
   aD( k) = 0;
   aDp(k) = 0;
 end

% coefficients aP and aPp

 aP  = VPC + (aU  + aD);
 aPp = VPC - (aUp + aDp);

% coefficient b

 b = dt * QS;

 switch upperBC_type
 case 'q'
   k = kU;
   b(k) = b(k) + (dtf *qs + dt1mf *qsn);
 end

 switch lowerBC_type
 case 'q'
   k = kL;
   b(k) = b(k) - (dtf *qb + dt1mf *qbn);
 end
