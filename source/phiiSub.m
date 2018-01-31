 function phi_i = phiiSub(Mtyp,Mw,phi_u)

% Finds the volume fraction of ice once Mw is known.

% Notation:

%   Mtyp    = material type                                 (M+1)
%   Mw      = mass of water (both phases) per unit volume   (M+1)
%   phi_u   = volume fraction of unfrozen water             (M+1)
%   phi_i   = volume fraction of ice                        (M+1)
% ______________________________________________

% Written by:

%   Gary Clow
%   Institute of Arctic and Alpine Research
%   University of Colorado
%   Boulder, Colorado USA
%   Email: gary.clow@colorado.edu
% ______________________________________________

% parameters

 rhoi = 917;            % density of ice
 rhou = 1000;           % density of water

% pre-allocate arrays

 phi_i = zeros(size(phi_u));

% find which material types are present

 [Mtypes,m] = which_Mtypes(Mtyp);

% > Find volume fraction of ice

 for j=1:m      % loop through unique material types

   L = Mtyp == Mtypes(j);

   switch fix(Mtypes(j))
   case 1                           % fixed properties (used for testing)
     phi_i(L) = 0;

   case {2,3}                       % pure ice
     phi_i(L) = 1;

   case {4,5,6,7,8,9,10,11,12,13,14,15,20,21,22,23,24,25}   % rocks and organic-rich materials

     phi_i(L)  = (Mw(L) - rhou *phi_u(L)) / rhoi;
     LL        = phi_i < 0;
     phi_i(LL) = 0;
   end
 end
