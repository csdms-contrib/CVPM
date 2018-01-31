 function [KeX,KeZ] = Keff_XZ(K,varepX,varepZ)

% Finds the effective conductivity at the CV interfaces for the XZ
% coordinate system.

% Notation:

%   K      = conductivity at CV grid points                 (M+1,N+1)
%   varepX = fractional distance horizontal direction       (N+1)
%   varepZ = fractional distance vertical direction         (M+1)
%   KeX    = effective conductivity across X-interfaces     (M+1,N+1)
%   KeZ    = effective conductivity across Z-interfaces     (M+1,N+1)
% ______________________________________________

% Written by:

%   Gary Clow
%   Institute of Arctic and Alpine Research
%   University of Colorado
%   Boulder, Colorado USA
%   Email: gary.clow@colorado.edu
% ______________________________________________

 N = length(varepX) - 1;
 M = length(varepZ) - 1;
 
% pre-allocate arrays

 KeX = NaN*ones(size(K));
 KeZ = NaN*ones(size(K));

% find effective conductivity at X-interfaces

 for k=1:M                      % do one row (k) at a time
   KX       = K(k,:);
   KeX(k,:) = Keff(KX,varepX);
 end
 KeX(M+1,:) = KeX(M,:);

% find effective conductivity at Z-interfaces

 for j=1:N                      % do one column (j) at a time
   KZ       = K(:,j);
   KeZ(:,j) = Keff(KZ,varepZ);
 end
 KeZ(:,N+1) = KeZ(:,N);
