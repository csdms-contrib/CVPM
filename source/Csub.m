 function [rhocp,C] = Csub(T,Mtyp,rhom,cpm0,phi,phi_i,phi_u,dphiudT)

% Finds the volumetric heat capacity C for each CV, including both latice
% vibration (rhocp) and latent-heat effects. For convenience, C is also 
% defined on the boundaries where it is set equal to that of the adjacent CVs.

% Notation:

%   T       = temperature (C)                       (M+1)
%   Mtyp    = material type                         (M+1)
%   rhom    = density of matrix material            (M+1)
%   cpm0    = matrix specific heat at 20 C          (M+1)
%   phi     = porosity                              (M+1)
%   phi_i   = volume fraction of ice                (M+1)
%   phi_u   = volume fraction of unfrozen water     (M+1)
%   dphiudT = d(phi_u)/dT                           (M+1)
%   C       = volumetric heat capacity              (M+1)
% ______________________________________________

% Written by:

%   Gary Clow
%   Institute of Arctic and Alpine Research
%   University of Colorado
%   Boulder, Colorado USA
%   Email: gary.clow@colorado.edu
% ______________________________________________

 M = length(T) - 1;

% pre-allocate arrays

 rhocp = NaN*ones(size(T));
 C     = NaN*ones(size(T));

% find which material types are present

 [Mtypes,m] = which_Mtypes(Mtyp);

% find rhocp and C for each material type present

 for j=1:m
   L = Mtyp == Mtypes(j);

   switch fix(Mtypes(j))
   case 1                           % fixed properties

     rhocp(L) = rhom(L) .* cpm0(L);
     C(L)     = rhocp(L);

   case 2                           % linearized ice

     rhocp(L) = C_linIce(T(L));
     C(L)     = rhocp(L);

   case 3                           % ice

     rhocp(L) = C_ice(T(L));
     C(L)     = rhocp(L);

   case {4,5,6,7,8,9,10,11,12,13,14,15}     % rocks and soils

     [rhocp(L),C(L)] = C_rock(T(L),Mtyp(L),rhom(L),cpm0(L),phi(L),phi_i(L),phi_u(L),dphiudT(L));

   case {20,21,22,23,24,25}         % organic-rich material

     [rhocp(L),C(L)] = C_veg(T(L),Mtyp(L),rhom(L),cpm0(L),phi(L),phi_i(L),phi_u(L),dphiudT(L));

   case {30,31,32,33,34,35}         % borehole fluids

     rhocp(L) = rhom(L) .* cpm0(L);
     C(L)     = rhocp(L);

   case {40,41,42,43,44}            % metals

     rhocp(L) = rhom(L) .* cpm0(L);
     C(L)     = rhocp(L);
   end
 end

% define C on the boundaries

 rhocp(1)   = rhocp(2);
 rhocp(M+1) = rhocp(M);
 C(1)       = C(2);
 C(M+1)     = C(M);
