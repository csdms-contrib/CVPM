 function [aU,aUp,aD,aDp,aW,aWp,aE,aEp,aP,aPp,b] ...
    = DEcoefs_XZ(dt,f,Dx,Dz,Ax,Az,VP,C,KeX,KeZ,QS,BCtypes,q,qn)

% Finds the discretization coefficients for the 2-D cylindrical CVPM model.

% Left X boundary condition is set by,

%   xleftBC_type:   'T' prescribed temperature
%                   'q' prescribed heat flux

% Right X boundary condition is set by,

%   xrighBC_type:	'T' prescribed temperature
%                   'q' prescribed heat flux

% Lower boundary condition is set by,

%   lowerBC_type:	'T' prescribed temperature
%                   'q' prescribed heat flux

% Notation:

%   dt                      time step                           (scalar)
%   f                       implicit/explicit factor            (scalar)
%   Dx                      horizontal CV width                 (N+1)
%   Dz                      vertical CV width                   (M+1)
%   Ax                      Dz/dx                               (M+1,N+1)
%   Az                      Dx/dz                               (M+1,N+1)
%   VP                      Dx*Dz                               (M+1,N+1)
%   C                       volumetric heat capacity            (M+1,N+1)
%   KeX                     effective conduct., X-interfaces    (M+1,N+1)                               (M+1,N+1)
%   KeZ                     effective conduct., Z-interfaces    (M+1,N+1)
%   QS                      spatially integrated source term    (M+1,N+1)
%   BCtypes                 inner/outer/upper/lower BC types    (strings)
%   q                       heat fluxes, this time step         (4)
%   qn                      heat fluxes, last time step         (4)
%   aW,aWp,aE,aEp,aU,aUp,   DE coefficients                     (M,N)
%   aD,aDp,aP,aPp,b          "      "                           (M,N)
%   j                       horizontal index
%   k                       vertical index

% Notes:

%   The source term is currently assumed to be independent of x.
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

% pre-allocate arrays

 aU  = zeros(M+1,N+1);
 aUp = zeros(M+1,N+1);
 aD  = zeros(M+1,N+1);
 aDp = zeros(M+1,N+1);
 aW  = zeros(M+1,N+1);
 aWp = zeros(M+1,N+1);
 aE  = zeros(M+1,N+1);
 aEp = zeros(M+1,N+1);
 b   = zeros(M+1,N+1);

% unpack BCs

 upperBC_type = BCtypes{1};
 lowerBC_type = BCtypes{2};
 xleftBC_type = BCtypes{3};
 xrighBC_type = BCtypes{4};

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

% indices of left/right/upper/lower CVs

 kU = 2;        % upper CV
 kL = M;        % lower CV
 jI = 2;        % left CV
 jO = N;        % right CV

% define the following for convenience

 g       = 4/3;
 dtf     = dt * f;
 dt1mf   = dt * (1-f);

 Dxdtf   = Dx * dtf;
 Dxdt1mf = Dx * dt1mf;
 Dzdtf   = Dz * dtf;
 Dzdt1mf = Dz * dt1mf;

 VPC     = VP .* C;
 AKeX    = Ax .* KeX;
 AKeZ    = Az .* KeZ;

% coefficients aU,aUp,aD,aDp

 for k=kU:kL
   aU( k,:) = dtf   * AKeZ(k  ,:);
   aUp(k,:) = dt1mf * AKeZ(k  ,:);
   aD( k,:) = dtf   * AKeZ(k+1,:);
   aDp(k,:) = dt1mf * AKeZ(k+1,:);
 end

 k = kU;
 switch upperBC_type
 case 'T'
   aU( k,:) = g * aU( k,:);
   aUp(k,:) = g * aUp(k,:);
   aD( k,:) = g * aD( k,:);
   aDp(k,:) = g * aDp(k,:);
 case 'q'
   aU( k,:) = 0;
   aUp(k,:) = 0;
 end

 k = kL;
 switch lowerBC_type
 case 'T'
   aU( k,:) = g * aU( k,:);
   aUp(k,:) = g * aUp(k,:);
   aD( k,:) = g * aD( k,:);
   aDp(k,:) = g * aDp(k,:);
 case 'q'
   aD( k,:) = 0;
   aDp(k,:) = 0;
 end

% coefficients aW,aWp,aE,aEp

 for j=jI:jO
   aW( :,j) = dtf   * AKeX(:,j);
   aWp(:,j) = dt1mf * AKeX(:,j);
   aE( :,j) = dtf   * AKeX(:,j+1);
   aEp(:,j) = dt1mf * AKeX(:,j+1);
 end

 j = jI;
 switch xleftBC_type
 case 'T'
   aW( :,j) = g * aW( :,j);
   aWp(:,j) = g * aWp(:,j);
   aE( :,j) = g * aE( :,j);
   aEp(:,j) = g * aEp(:,j);
 case 'q'
   aW( :,j) = 0;
   aWp(:,j) = 0;
 end

 j = jO;
 switch xrighBC_type
 case 'T'
   aW( :,j) = g * aW( :,j);
   aWp(:,j) = g * aWp(:,j);
   aE( :,j) = g * aE( :,j);
   aEp(:,j) = g * aEp(:,j);
 case 'q'
   aE( :,j) = 0;
   aEp(:,j) = 0;
 end

% coefficients aP and aPp

 aP  = VPC + (aU  + aD  + aW  + aE);
 aPp = VPC - (aUp + aDp + aWp + aEp);

% coefficient b

 b = dt * repmat(Dx,[M+1 1]) .* QS;

 switch upperBC_type
 case 'q'
   k = kU;
   b(k,:) = b(k,:) + (Dxdtf .*qs + Dxdt1mf .*qsn);
 end
 switch lowerBC_type
 case 'q'
   k = kL;
   b(k,:) = b(k,:) - (Dxdtf .*qb + Dxdt1mf .*qbn);
 end

 switch xleftBC_type
 case 'q'
   j = jI;
   b(:,j) = b(:,j) + (Dzdtf .*qa + Dzdt1mf .*qan);
 end
 switch xrighBC_type
 case 'q'
   j = jO;
   b(:,j) = b(:,j) - (Dzdtf .*qo + Dzdt1mf .*qon);
 end

% discard unused portion of arrays

 aU  = aU( 1:M,1:N);
 aUp = aUp(1:M,1:N);
 aD  = aD( 1:M,1:N);
 aDp = aDp(1:M,1:N);
 aW  = aW( 1:M,1:N);
 aWp = aWp(1:M,1:N);
 aE  = aE( 1:M,1:N);
 aEp = aEp(1:M,1:N);
 aP  = aP( 1:M,1:N);
 aPp = aPp(1:M,1:N);
 b   = b(  1:M,1:N);
