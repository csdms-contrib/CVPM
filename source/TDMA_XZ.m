 function T = TDMA_XZ(A,B,C,D,E,Ap,Bp,Cp,Dp,Ep,b,Tn,BCtypes,Tbc)

% TDMA algorithm for the 2-D cartesian CVPM model.
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

% This version employs TDMA along vertical columns, sweeping from top to bottom.

% Currently allowed boundary conditions:

%   upperBC_type:	'T'     prescribed temperature
%   lowerBC_type:	'T'     prescribed temperature
%                   'q'     prescribed heat flux
%   xleftBC_type:	'T'     prescribed temperature
%                   'q'     prescribed heat flux
%   xrighBC_type:	'T'     prescribed temperature
%                   'q'     prescribed heat flux
% ______________________________________________

 M = size(Tn,1) - 1;
 N = size(Tn,2) - 1;

% pre-allocate arrays

 T = zeros(size(Tn));

% unpack BCs

 lowerBC_type = BCtypes{2};
 xleftBC_type = BCtypes{3};
 xrighBC_type = BCtypes{4};

 Ts = Tbc{1};
 Tb = Tbc{2};
 Ta = Tbc{3};
 To = Tbc{4};

% temperatures on upper boundary

 T(1,:) = Ts;

% > Perform TDMA for each CV column

 for j=2:N

   F = zeros(M,1);
   P = zeros(M+1,1);
   Q = zeros(M+1,1);

% coefficient F

   F(2:M) = Ap(2:M,j).*Tn(1:M-1,j) + Cp(2:M,j).*Tn(3:M+1,j) ...
          + Dp(2:M,j).*Tn(2:M,j-1) + Ep(2:M,j).*Tn(2:M,j+1) ...
          + Bp(2:M,j).*Tn(2:M,j)   + b(2:M,j);

   switch lowerBC_type
   case 'q'
     F(M) = Ap(M,j)*Tn(M-1,j) + Dp(M,j)*Tn(M,j-1) + Ep(M,j)*Tn(M,j+1) ...
          + Bp(M,j)*Tn(M,j)   + b(M,j);
   end

% recursion coefficients P,Q

   switch lowerBC_type
   case 'T'
     P(M+1) = 0;
     Q(M+1) = Tb(j);
   case 'q'
     P(M+1) = 0;
     Q(M+1) = 0;
   end

   for k=M:-1:2
     fac  =  B(k,j) - C(k,j)*P(k+1);
     P(k) =  A(k,j) / fac;
     Q(k) = (C(k,j)*Q(k+1) + D(k,j)*T(k,j-1) + E(k,j)*Tn(k,j+1) + F(k)) / fac;
   end

% temperature at internal grid points

   for k=2:M
     T(k,j) = P(k)*T(k-1,j) + Q(k);
   end
 end

% > Set temperatures on other boundaries as appropriate

% lower boundary

 switch lowerBC_type
 case 'T'
   T(M+1,:) = Tb;
 end

% x-left boundary

 switch xleftBC_type
 case 'T'
   T(:,1) = Ta;
 end

% x-right boundary

 switch xrighBC_type
 case 'T'
   T(:,N+1) = To;
 end
