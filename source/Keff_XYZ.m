 function [KeX,KeY,KeZ] = Keff_XYZ(K,varepX,varepY,varepZ)

% Finds the effective conductivity at the CV interfaces for the XYZ
% coordinate system.

% Notation:

%   K      = conductivity at CV grid points                 (M+1,N+1,L+1)
%   varepX = fractional distance horizontal direction       (N+1)
%   varepZ = fractional distance vertical direction         (M+1)
%   KeX    = effective conductivity across X-interfaces     (M+1,N+1,L+1)
%   KeY    = effective conductivity across Y-interfaces     (M+1,N+1,L+1)
%   KeZ    = effective conductivity across Z-interfaces     (M+1,N+1,L+1)
% ______________________________________________

% Written by:

%   Gary Clow
%   Institute of Arctic and Alpine Research
%   University of Colorado
%   Boulder, Colorado USA
%   Email: gary.clow@colorado.edu
% ______________________________________________

 N = length(varepX) - 1;
 L = length(varepY) - 1;
 M = length(varepZ) - 1;
 
% pre-allocate arrays

 KeX = NaN*ones(size(K));
 KeY = NaN*ones(size(K));
 KeZ = NaN*ones(size(K));

% find effective conductivity at X-interfaces

 for jy=1:L
   for k=1:M                      % do one row (k) at a time
     KX          = K(k,:,jy);
     KeX(k,:,jy) = Keff(KX,varepX);
   end
 end
 KeX(M+1,:,:) = KeX(M,:,:);
 KeX(:,:,L+1) = KeX(:,:,L);

% find effective conductivity at Y-interfaces

 for j=1:N
   for k=1:M                      % do one row (k) at a time
     KY         = K(k,j,:);
     KeY(k,j,:) = Keff(KY,varepY);
   end
 end
 KeY(M+1,:,:) = KeY(M,:,:);
 KeY(:,N+1,:) = KeY(:,N,:);

% find effective conductivity at Z-interfaces

 for jy=1:L
   for j=1:N                      % do one column (j) at a time
     KZ          = K(:,j,jy);
     KeZ(:,j,jy) = Keff(KZ,varepZ);
   end
 end
 KeZ(:,N+1,:) = KeZ(:,N,:);
 KeZ(:,:,L+1) = KeZ(:,:,L);
