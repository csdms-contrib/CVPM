 function phi_u = phiuSub(Mtyp,phi_w,lambda,r,DeltaT)

% Finds the volume fraction of unfrozen water during the initial
% setup by CPS.

% Notation:

%   T       = temperature (C)                               (M+1)
%   Mtyp    = material type                                 (M+1)
%   rhom    = matrix conductivity                           (M+1)
%   phi     = porosity                                      (M+1)
%   phi_w   = volume fraction of water (solid & liquid)     (M+1)
%   lambda  = interfacial melting parameter                 (M+1)
%   r       = effective radius of matrix particles          (M+1)
%   DeltaT  = Tf - T                                        (M+1)
%   phi_u   = volume fraction of unfrozen water             (M+1)
% ______________________________________________

% Written by:

%   Gary Clow
%   Institute of Arctic and Alpine Research
%   University of Colorado
%   Boulder, Colorado USA
%   Email: gary.clow@colorado.edu
% ______________________________________________

% parameters

 PC_opt = 'fcc';

% pre-allocate arrays

 phi_u = zeros(size(phi_w));

% find which material types are present

 [Mtypes,m] = which_Mtypes(Mtyp);

% > Find volume fraction of unfrozen water

 for j=1:m      % loop through unique material types

   L = Mtyp == Mtypes(j);

   switch fix(Mtypes(j))
   case 1                           % fixed properties (used for testing)
     phi_u(L) = 0;

   case {2,3}                       % pure ice
     phi_u(L) = 0;

   case {4,5,6,7,8,9,10,11,12,13,14,15,20,21,22,23,24,25}   % rocks and organic-rich materials

     phi_u(L) = phiu_pore(phi_w(L),lambda(L),r(L),DeltaT(L),PC_opt);
   end
 end
