 function T = TDMA_RZh(A,B,C,D,E,Ap,Bp,Cp,Dp,Ep,b,Tn,BCtypes,Tbc)

% TDMA algorithm for the 2-D cylindrical CVPM model.
% This version employs TDMA along horizontal rows, sweeping from left to right.

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

% > Perform TDMA for each CV row

 for k=2:M

   F = zeros(1,N);
   P = zeros(1,N+1);
   Q = zeros(1,N+1);

% coefficient F

   F(2:N) = Ap(k,2:N).*Tn(k-1,2:N) + Bp(k,2:N).*Tn(k,2:N) ...
          + Cp(k,2:N).*Tn(k+1,2:N) + Dp(k,2:N).*Tn(k,1:N-1) ...
                                   + Ep(k,2:N).*Tn(k,3:N+1) + b(k,2:N);

   switch innerBC_type
   case 'none'
     F(1) = Ap(k,1)*Tn(k-1,1) + Bp(k,1)*Tn(k,1) + Cp(k,1)*Tn(k+1,1) ...
          + Ep(k,1)*Tn(k,2) + b(k,1);
   end

   switch outerBC_type
   case 'q'
     F(N) = Ap(k,N)*Tn(k-1,N) + Bp(k,N)*Tn(k,N) + Cp(k,N)*Tn(k+1,N) ...
          + Dp(k,N)*Tn(k,N-1) + b(k,N);
   end

% recursion coefficients P,Q

   switch outerBC_type
   case 'T'
     P(N+1) = 0;
     Q(N+1) = To(k);
   case 'q'
     P(N+1) = 0;
     Q(N+1) = 0;
   end

   for j=N:-1:2
     fac  =  B(k,j) - E(k,j)*P(j+1);
     P(j) =  D(k,j) / fac;
     Q(j) = (A(k,j)*T(k-1,j) + C(k,j)*Tn(k+1,j) + E(k,j)*Q(j+1) + F(j)) / fac;
   end

% inner BC

   switch innerBC_type
   case 'T'
     T(k,1) = Ta(k);
   case 'none'
     fac  =  B(k,1) - E(k,1)*P(2);
     P(1) =  0;
     Q(1) = (A(k,1)*T(k-1,1) + C(k,1)*Tn(k+1,1) + E(k,1)*Q(2) + F(1)) / fac;
     T(k,1) = Q(1);                    % this grid point is at the very center
   end

% temperature at internal grid points

   for j=2:N
     T(k,j) = P(j)*T(k,j-1) + Q(j);
   end
 end

% > Set temperatures on other boundaries as appropriate

% lower boundary

 switch lowerBC_type
 case 'T'
   T(M+1,:) = Tb;
   T(M+1,1) = Ta(1);	% assign wall temperature to corner point
 end

% outer boundary

 switch outerBC_type
 case 'T'
   T(:,N+1) = To;
 end
