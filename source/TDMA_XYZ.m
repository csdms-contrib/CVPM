 function T = TDMA_XYZ(A,B,C,D,E,G,H,Ap,Bp,Cp,Dp,Ep,Gp,Hp,b,Tn,BCtypes,Tbc)

% TDMA algorithm for the 3-D cartesian CVPM model.
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

% This version employs TDMA along vertical columns.

% Currently allowed boundary conditions:

% Upper BC:   prescribed temperature Ts(x,y,t) at earth's surface
% Lower BC:   prescribed heat flux   qb(x,y,t) at lower boundary
% x Left  BC: prescribed temperature Ta(y,z,t) or flux at x-left boundary
% x Right BC: prescribed temperature To(y,z,t) or flux at x-right boundary
% y Left  BC: prescribed temperature Ta(y,z,t) or flux at y-left boundary
% y Right BC: prescribed temperature To(y,z,t) or flux at y-right boundary
% ______________________________________________

 M = size(Tn,1) - 1;
 N = size(Tn,2) - 1;
 L = size(Tn,3) - 1;

% pre-allocate arrays

 T = zeros(size(Tn));

% unpack BCs

 upperBC_type = BCtypes{1};
 lowerBC_type = BCtypes{2};
 xleftBC_type = BCtypes{3};
 xrighBC_type = BCtypes{4};
 yleftBC_type = BCtypes{5};
 yrighBC_type = BCtypes{6};

 Ts = Tbc{1};
 Tb = Tbc{2};
 Ta = Tbc{3};
 To = Tbc{4};
 Tc = Tbc{5};
 Td = Tbc{6};

% temperatures prescribed on upper boundary

 T(1,:,:) = Ts;

% > Perform TDMA for each CV column

 for i=2:L
   for j=2:N
   
     F = zeros(M,  1);
     P = zeros(M+1,1);
     Q = zeros(M+1,1);

%   coefficient F

     F(2:M) = Ap(2:M,j,i).*Tn(1:M-1,j,i) + Cp(2:M,j,i).*Tn(3:M+1,j,i) ...
            + Dp(2:M,j,i).*Tn(2:M,j-1,i) + Ep(2:M,j,i).*Tn(2:M,j+1,i) ...
            + Gp(2:M,j,i).*Tn(2:M,j,i-1) + Hp(2:M,j,i).*Tn(2:M,j,i+1) ...
            + Bp(2:M,j,i).*Tn(2:M,j,i)   +  b(2:M,j,i);

     switch lowerBC_type
     case 'q'
       F(M) = Ap(M,j,i)*Tn(M-1,j,i) ...
            + Dp(M,j,i)*Tn(M,j-1,i) + Ep(M,j,i)*Tn(M,j+1,i) ...
            + Gp(M,j,i)*Tn(M,j,i-1) + Hp(M,j,i)*Tn(M,j,i+1) ...
            + Bp(M,j,i)*Tn(M,j,i)   +  b(M,j,i);
     end

%   recursion coefficients P,Q

     switch lowerBC_type
     case 'T'
       P(M+1) = 0;
       Q(M+1) = Tb(j,i);
     case 'q'
       P(M+1) = 0;
       Q(M+1) = 0;
     end

     for k=M:-1:2
       fac  =  B(k,j,i) - C(k,j,i)*P(k+1);
       P(k) =  A(k,j,i) / fac;
       Q(k) = (C(k,j,i)*Q(k+1) + D(k,j,i)*T(k,j-1,i) + E(k,j,i)*Tn(k,j+1,i) ...
             + G(k,j,i)*T(k,j,i-1) + H(k,j,i)*Tn(k,j,i+1) + F(k)) / fac;
     end

%   temperature at internal grid points

     for k=2:M
       T(k,j,i) = P(k)*T(k-1,j,i) + Q(k);
     end
   end
 end

% temperature at lower boundary

 switch lowerBC_type
 case 'T'
   T(M+1,:,:) = Tb;
 end

% temperature at x-left boundary

 switch xleftBC_type
 case 'T'
   T(:,1,:) = Ta;
 end

% temperature at x-right boundary

 switch xrighBC_type
 case 'T'
   T(:,N+1,:) = To;
 end

% temperature at y-left boundary

 switch yleftBC_type
 case 'T'
   T(:,:,1) = Tc;
 end

% temperature at y-right boundary

 switch yrighBC_type
 case 'T'
   T(:,:,L+1) = Td;
 end
