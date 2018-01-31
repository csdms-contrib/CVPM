 function [aU,aUp,aD,aDp,aW,aWp,aE,aEp,aS,aSp,aN,aNp,aP,aPp,b] ...
      = DEcoefs_XYZ(dt,f,Dx,Dy,Dz,Ax,Ay,Az,VP,C,KeX,KeY,KeZ,QS,BCtypes,q,qn)

% Finds the discretization coefficients for the 3-D cylindrical CVPM model.

% Left X boundary condition is set by,

%   xleftBC_type:   'T' prescribed temperature
%                   'q' prescribed heat flux

% Right X boundary condition is set by,

%   xrighBC_type:	'T' prescribed temperature
%                   'q' prescribed heat flux

% Left Y boundary condition is set by,

%   yleftBC_type:	'T' prescribed temperature
%                   'q' prescribed heat flux

% Right Y boundary condition is set by,

%   yrighBC_type:	'T' prescribed temperature
%                   'q' prescribed heat flux

% Lower boundary condition is set by,

%   lowerBC_type:	'T' prescribed temperature
%                   'q' prescribed heat flux

% Notation:

%   dt                      time step                           (scalar)
%   f                       implicit/explicit factor            (scalar)
%   Dx                      horizontal CV width                 (N+1)
%   Dy                      horizontal CV width                 (L+1)
%   Dz                      vertical CV width                   (M+1)
%   Ax                      Dy*Dz/dx                            (M+1,N+1,L+1)
%   Ay                      Dx*Dz/dy                            (M+1,N+1,L+1)
%   Az                      Dx*Dy/dz                            (M+1,N+1,L+1)
%   VP                      Dx*Dy*Dz                            (M+1,N+1,L+1)
%   C                       volumetric heat capacity            (M+1,N+1,L+1)
%   KeX                     effective conduct., X-interfaces    (M+1,N+1,L+1)
%   KeY                     effective conduct., Y-interfaces    (M+1,N+1,L+1)
%   KeZ                     effective conduct., Z-interfaces    (M+1,N+1,L+1)
%   QS                      spatially integrated source term    (M+1,N+1,L+1)
%   BCtypes                 inner/outer/upper/lower BC types    (strings)
%   q                       heat fluxes, this time step         (6)
%   qn                      heat fluxes, last time step         (6)
%   aW,aWp,aE,aEp,aS,aSp,   DE coefficients                     (M,N,L)
%   aN,aNp,aU,aUp,aD,aDp,aP,aPp,b                               (M,N,L)
%   j                       x horizontal index
%   jy                      y horizontal index
%   k                       vertical index

% Notes:

%   The source term is assumed to be independent of x and y.
% ______________________________________________

% Written by:

%   Gary Clow
%   Institute of Arctic and Alpine Research
%   University of Colorado
%   Boulder, Colorado USA
%   Email: gary.clow@colorado.edu
% ______________________________________________

 M = size(C,1) - 1;
 N = size(C,2) - 1;
 L = size(C,3) - 1;

% pre-allocate arrays

 aU  = zeros(M+1,N+1,L+1);
 aUp = zeros(M+1,N+1,L+1);
 aD  = zeros(M+1,N+1,L+1);
 aDp = zeros(M+1,N+1,L+1);
 aW  = zeros(M+1,N+1,L+1);
 aWp = zeros(M+1,N+1,L+1);
 aE  = zeros(M+1,N+1,L+1);
 aEp = zeros(M+1,N+1,L+1);
 aS  = zeros(M+1,N+1,L+1);
 aSp = zeros(M+1,N+1,L+1);
 aN  = zeros(M+1,N+1,L+1);
 aNp = zeros(M+1,N+1,L+1);
 b   = zeros(M+1,N+1,L+1);

% unpack BCs

 upperBC_type = BCtypes{1};
 lowerBC_type = BCtypes{2};
 xleftBC_type = BCtypes{3};
 xrighBC_type = BCtypes{4};
 yleftBC_type = BCtypes{5};
 yrighBC_type = BCtypes{6};

 switch upperBC_type
 case 'q'
   qs  = q{ 1};
   qsn = qn{1};
 end
 switch lowerBC_type
 case 'q'
   qb  = q{ 2};
   qbn = qn{2};
 end

 switch xleftBC_type
 case 'q'
   qa  = q{ 3};
   qan = qn{3};
 end
 switch xrighBC_type
 case 'q'
   qo  = q{ 4};
   qon = qn{4};
 end

 switch yleftBC_type
 case 'q'
   qc  = q{ 5};
   qcn = qn{5};
 end
 switch yrighBC_type
 case 'q'
   qd  = q{ 6};
   qdn = qn{6};
 end

% indices of left/right/upper/lower CVs

 kU  = 2;       % upper CV
 kL  = M;       % lower CV
 jI  = 2;       % x left CV
 jO  = N;       % x right CV
 jyI = 2;       % y left CV
 jyO = L;       % y right CV

% define the following for convenience

 g         = 4/3;
 dtf       = dt * f;
 dt1mf     = dt * (1-f);

 DxDz      = repmat(Dx, M+1,1) .* repmat(Dz,1,N+1);
 DyDz      = repmat(Dy, M+1,1) .* repmat(Dz,1,L+1);
 DxDy      = repmat(Dx',1,L+1) .* repmat(Dy,N+1,1);

 DxDzdtf   = DxDz * dtf;
 DxDzdt1mf = DxDz * dt1mf;
 DyDzdtf   = DyDz * dtf;
 DyDzdt1mf = DyDz * dt1mf;
 DxDydtf   = DxDy * dtf;
 DxDydt1mf = DxDy * dt1mf;

 VPC     = VP .* C;
 AKeX    = Ax .* KeX;
 AKeY    = Ay .* KeY;
 AKeZ    = Az .* KeZ;

% coefficients aU,aUp,aD,aDp

 for k=kU:kL
   aU( k,:,:) = dtf   * AKeZ(k  ,:,:);
   aUp(k,:,:) = dt1mf * AKeZ(k  ,:,:);
   aD( k,:,:) = dtf   * AKeZ(k+1,:,:);
   aDp(k,:,:) = dt1mf * AKeZ(k+1,:,:);
 end

 k = kU;
 switch upperBC_type
 case 'T'
   aU( k,:,:) = g * aU( k,:,:);
   aUp(k,:,:) = g * aUp(k,:,:);
   aD( k,:,:) = g * aD( k,:,:);
   aDp(k,:,:) = g * aDp(k,:,:);
 case 'q'
   aU( k,:,:) = 0;
   aUp(k,:,:) = 0;
 end

 k = kL;
 switch lowerBC_type
 case 'T'
   aU( k,:,:) = g * aU( k,:,:);
   aUp(k,:,:) = g * aUp(k,:,:);
   aD( k,:,:) = g * aD( k,:,:);
   aDp(k,:,:) = g * aDp(k,:,:);
 case 'q'
   aD( k,:,:) = 0;
   aDp(k,:,:) = 0;
 end

% coefficients aW,aWp,aE,aEp

 for j=jI:jO
   aW( :,j,:) = dtf   * AKeX(:,j,  :);
   aWp(:,j,:) = dt1mf * AKeX(:,j,  :);
   aE( :,j,:) = dtf   * AKeX(:,j+1,:);
   aEp(:,j,:) = dt1mf * AKeX(:,j+1,:);
 end

 j = jI;
 switch xleftBC_type
 case 'T'
   aW( :,j,:) = g * aW( :,j,:);
   aWp(:,j,:) = g * aWp(:,j,:);
   aE( :,j,:) = g * aE( :,j,:);
   aEp(:,j,:) = g * aEp(:,j,:);
 case 'q'
   aW( :,j,:) = 0;
   aWp(:,j,:) = 0;
 end

 j = jO;
 switch xrighBC_type
 case 'T'
   aW( :,j,:) = g * aW( :,j,:);
   aWp(:,j,:) = g * aWp(:,j,:);
   aE( :,j,:) = g * aE( :,j,:);
   aEp(:,j,:) = g * aEp(:,j,:);
 case 'q'
   aE( :,j,:) = 0;
   aEp(:,j,:) = 0;
 end

% coefficients aS,aSp,aN,aNp

 for jy=jyI:jyO
   aS( :,:,jy) = dtf   * AKeY(:,:,jy);
   aSp(:,:,jy) = dt1mf * AKeY(:,:,jy);
   aN( :,:,jy) = dtf   * AKeY(:,:,jy+1);
   aNp(:,:,jy) = dt1mf * AKeY(:,:,jy+1);
 end

 jy = jyI;
 switch yleftBC_type
 case 'T'
   aS( :,:,jy) = g * aS( :,:,jy);
   aSp(:,:,jy) = g * aSp(:,:,jy);
   aN( :,:,jy) = g * aN( :,:,jy);
   aNp(:,:,jy) = g * aNp(:,:,jy);
 case 'q'
   aS( :,:,jy) = 0;
   aSp(:,:,jy) = 0;
 end

 jy = jyO;
 switch yrighBC_type
 case 'T'
   aS( :,:,jy) = g * aS( :,:,jy);
   aSp(:,:,jy) = g * aSp(:,:,jy);
   aN( :,:,jy) = g * aN( :,:,jy);
   aNp(:,:,jy) = g * aNp(:,:,jy);
 case 'q'
   aN( :,:,jy) = 0;
   aNp(:,:,jy) = 0;
 end

% coefficients aP and aPp

 aP  = VPC + (aU  + aD  + aW  + aE  + aS  + aN);
 aPp = VPC - (aUp + aDp + aWp + aEp + aSp + aNp);

% coefficient b

 Dyy = permute(Dy,[1 3 2]);
 b   = dt * repmat(Dx,[M+1 1 L+1]) .* repmat(Dyy,[M+1 N+1 1]) .* QS;

 switch upperBC_type
 case 'q'
   k = kU;
   A = DxDydtf .*qs + DxDydt1mf .*qsn;
   b(k,:,:) = b(k,:,:) + permute(A,[3 1 2]);
 end
 switch lowerBC_type
 case 'q'
   k = kL;
   A = DxDydtf .*qb + DxDydt1mf .*qbn;
   b(k,:,:) = b(k,:,:) - permute(A,[3 1 2]);
 end

 switch xleftBC_type
 case 'q'
   j = jI;
   A = DyDzdtf .*qa + DyDzdt1mf .*qan;
   b(:,j,:) = b(:,j,:) + permute(A,[1 3 2]);
 end
 switch xrighBC_type
 case 'q'
   j = jO;
   A = DyDzdtf .*qo + DyDzdt1mf .*qon;
   b(:,j,:) = b(:,j,:) - permute(A,[1 3 2]);
 end

 switch yleftBC_type
 case 'q'
   jy = jyI;
   A  = DxDzdtf .*qc + DxDzdt1mf .*qcn;
   b(:,:,jy) = b(:,:,jy) + permute(A,[1 2 3]);
 end
 switch yrighBC_type
 case 'q'
   jy = jyO;
   A  = DxDzdtf .*qd + DxDzdt1mf .*qdn;
   b(:,:,jy) = b(:,:,jy) + permute(A,[1 2 3]);
 end

% discard unused portion of arrays

 aU  = aU( 1:M,1:N,1:L);
 aUp = aUp(1:M,1:N,1:L);
 aD  = aD( 1:M,1:N,1:L);
 aDp = aDp(1:M,1:N,1:L);
 aW  = aW( 1:M,1:N,1:L);
 aWp = aWp(1:M,1:N,1:L);
 aE  = aE( 1:M,1:N,1:L);
 aEp = aEp(1:M,1:N,1:L);
 aS  = aS( 1:M,1:N,1:L);
 aSp = aSp(1:M,1:N,1:L);
 aN  = aN( 1:M,1:N,1:L);
 aNp = aNp(1:M,1:N,1:L);
 aP  = aP( 1:M,1:N,1:L);
 aPp = aPp(1:M,1:N,1:L);
 b   = b(  1:M,1:N,1:L);
