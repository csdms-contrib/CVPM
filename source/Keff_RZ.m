 function [KeR,KeZ] = Keff_RZ(K,varepR,varepZ)

% Finds the effective conductivity at the CV interfaces for the RZ
% coordinate system.

% Notation:

%   K      = conductivity at CV grid points                 (M+1,N+1)
%   varepR = fractional distance radial direction           (N+1)
%   varepZ = fractional distance vertical direction         (M+1)
%   KeR    = effective conductivity across R-interfaces     (M+1,N+1)
%   KeZ    = effective conductivity across Z-interfaces     (M+1,N+1)
% ______________________________________________

% Written by:

%   Gary Clow
%   Institute of Arctic and Alpine Research
%   University of Colorado
%   Boulder, Colorado USA
%   Email: gary.clow@colorado.edu
% ______________________________________________

 N = length(varepR) - 1;
 M = length(varepZ) - 1;
 
% pre-allocate array

 KeR = NaN*ones(size(K));
 KeZ = NaN*ones(size(K));

% find effective conductivity at R-interfaces

 for k=1:M                      % do one row (k) at a time
   KR       = K(k,:);
   KeR(k,:) = Keff(KR,varepR);
 end
 KeR(M+1,:) = KeR(M,:);

% find effective conductivity at Z-interfaces

 for j=1:N                      % do one column (j) at a time
   KZ       = K(:,j);
   KeZ(:,j) = Keff(KZ,varepZ);
 end
 KeZ(:,N+1) = KeZ(:,N);
