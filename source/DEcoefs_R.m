 function [aW,aWp,aE,aEp,aP,aPp,b] ...
    = DEcoefs_R(dt,f,rf,Ar,VP,C,Ke,QS,BCtypes,q,qn)

% Finds the discretization coefficients for 1-D radial CVPM model.

% Inner boundary condition is set by,

%   innerBC_type:   'T'     prescribed temperature
%                   'q'     prescribed heat flux
%                   'none'  there is no inner boundary

% Outer boundary condition is set by,

%   outerBC_type:	'T' prescribed temperature
%                   'q' prescribed heat flux

% Notation:

%   dt                      time step                           (scalar)
%   f                       implicit/explicit factor            (scalar)
%   rf                      CV interfaces                       (N+1)
%   Ar                      rf/dr                               (N+1)
%   VP                      lambda                              (N+1)
%   C                       volumetric heat capacity            (N+1)
%   Ke                      effective conduct., R-interfaces    (N+1)
%   QS                      spatially integrated source term    (N+1)
%   BCtypes                 inner/outer BC types                (strings)
%   q                       heat fluxes, this time step         (2)
%   qn                      heat fluxes, last time step         (2)
%   aW,aWp,aE,aEp,aP,aPp,b  DE coefficients                     (N+1)
% ______________________________________________

% Written by:

%   Gary Clow
%   Institute of Arctic and Alpine Research
%   University of Colorado
%   Boulder, Colorado USA
%   Email: gary.clow@colorado.edu
% ______________________________________________

 N = length(C) - 1;

% pre-allocate arrays

 aW  = zeros(1,N+1);
 aWp = zeros(1,N+1);
 aE  = zeros(1,N+1);
 aEp = zeros(1,N+1);

% unpack BCs

 innerBC_type = BCtypes{1};
 outerBC_type = BCtypes{2};

 switch innerBC_type
 case 'q'
   qa  = q(1);
   qan = qn(1);
 end
 switch outerBC_type
 case 'q'
   qo  = q(2);
   qon = qn(2);
 end

% indices of inner and outer CVs

 switch innerBC_type
 case {'T','q'}
   jI = 2;      % inner CV
 case 'none'
   jI = 1;
 end
 jO = N;        % outer CV

% define the following for convenience

 g     = 4/3;
 dtf   = dt * f;
 dt1mf = dt * (1-f);
 VPC   = VP .* C;
 AKeR  = Ar .* Ke;

% coefficients aW,aWp,aE,aEp

 aW( jI:jO) = dtf   * AKeR(jI:jO);
 aWp(jI:jO) = dt1mf * AKeR(jI:jO);
 aE( jI:jO) = dtf   * AKeR(jI+1:jO+1);
 aEp(jI:jO) = dt1mf * AKeR(jI+1:jO+1);

 j = jI;
 switch innerBC_type
 case 'T'
   aW( j) = g * aW( j);
   aWp(j) = g * aWp(j);
   aE( j) = g * aE( j);
   aEp(j) = g * aEp(j);
 case {'q','none'}
   aW( j) = 0;
   aWp(j) = 0;
 end

 j = jO;
 switch outerBC_type
 case 'T'
   aW( j) = g * aW( j);
   aWp(j) = g * aWp(j);
   aE( j) = g * aE( j);
   aEp(j) = g * aEp(j);
 case 'q'
   aE( j) = 0;
   aEp(j) = 0;
 end

% coefficients aP and aPp

 aP  = VPC + (aW  + aE);
 aPp = VPC - (aWp + aEp);

% coefficient b

 b = dt * QS;

 switch innerBC_type
 case 'q'
   j = jI;
   b(j) = b(j) + (dtf *rf(j)*qa + dt1mf *rf(j)*qan);
 end

 switch outerBC_type
 case 'q'
   j = jO;
   b(j) = b(j) - (dtf *rf(j+1)*qo + dt1mf *rf(j+1)*qon);
 end
