 function T = TDMA_RZv(A,B,C,D,E,Ap,Bp,Cp,Dp,Ep,b,Tn,BCtypes,Tbc)

% TDMA algorithm for the 2-D cylindrical CVPM model.
% This version employs TDMA along vertical columns, sweeping from top to bottom.

% Currently allowed boundary conditions:

%   upperBC_type:	'T'     prescribed temperature
%   lowerBC_type:	'T'     prescribed temperature
%                   'q'     prescribed heat flux
%   innerBC_type:   'T'     prescribed temperature
%                   'none'  no boundary condition
%   outerBC_type:	'T'     prescribed temperature
%                   'q'     prescribed heat flux
% ______________________________________________

% Written by:

%   Gary Clow
%   Institute of Arctic and Alpine Research
%   University of Colorado
%   Boulder, Colorado USA
%   Email: gary.clow@colorado.edu
% ______________________________________________

 M = size(Tn,1) - 1;
 N = size(Tn,2) - 1;

% pre-allocate arrays

 T = zeros(size(Tn));

% unpack BCs

 lowerBC_type = BCtypes{2};
 innerBC_type = BCtypes{3};
 outerBC_type = BCtypes{4};

 Ts = Tbc{1};
 Tb = Tbc{2};
 Ta = Tbc{3};
 To = Tbc{4};

% temperatures on upper boundary

 T(1,:) = Ts;

 switch innerBC_type
 case 'T'
   T(1,1) = Ta(1);    % assign wall temperature to corner point
 end

% > Perform TDMA for each CV column

 switch innerBC_type
 case {'T','q'}
   jI = 2;
 case 'none'
   jI = 1;
 end

 for j=jI:N

   F = zeros(M,1);
   P = zeros(M+1,1);
   Q = zeros(M+1,1);

% coefficient F

   switch j
   case 1
     F(2:M) = Ap(2:M,j).*Tn(1:M-1,j) + Cp(2:M,j).*Tn(3:M+1,j) ...
            + Ep(2:M,j).*Tn(2:M,j+1) ...
            + Bp(2:M,j).*Tn(2:M,j)   + b(2:M,j);
   otherwise
     F(2:M) = Ap(2:M,j).*Tn(1:M-1,j) + Cp(2:M,j).*Tn(3:M+1,j) ...
            + Dp(2:M,j).*Tn(2:M,j-1) + Ep(2:M,j).*Tn(2:M,j+1) ...
            + Bp(2:M,j).*Tn(2:M,j)   + b(2:M,j);
   end

   switch lowerBC_type
   case 'q'
     switch j
     case 1
       F(M) = Ap(M,j)*Tn(M-1,j) + Ep(M,j)*Tn(M,j+1) ...
            + Bp(M,j)*Tn(M,j)   + b(M,j);
     otherwise
       F(M) = Ap(M,j)*Tn(M-1,j) + Dp(M,j)*Tn(M,j-1) + Ep(M,j)*Tn(M,j+1) ...
            + Bp(M,j)*Tn(M,j)   + b(M,j);
     end
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
     switch j
     case 1
       Q(k) = (C(k,j)*Q(k+1) + E(k,j)*Tn(k,j+1) + F(k)) / fac;
     otherwise
       Q(k) = (C(k,j)*Q(k+1) + D(k,j)*T(k,j-1) + E(k,j)*Tn(k,j+1) + F(k)) / fac;
     end
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

% inner boundary

 switch innerBC_type
 case 'T'
   T(:,1) = Ta;
   T(M+1,1) = Ta(1);	% assign wall temperature to corner point
 end

% outer boundary

 switch outerBC_type
 case 'T'
   T(:,N+1) = To;
 end
