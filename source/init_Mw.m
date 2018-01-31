 function Mw = init_Mw(Mtyp,phi_i,phi_u)

% Initializes the total water mass per unit volume.

% Notation:

%   Mtyp  = material type                           (M+1)
%   phi_i = volume fraction of ice                  (M+1)
%   phi_u = volume fraction of unfrozen water       (M+1)
%   Mw    = volumetric water mass                   (M+1)
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
 rhou = 1000;           % density of unfrozen water

% pre-allocate the array

 Mw = zeros(size(Mtyp));

% reset the water mass for rocks, soils, and organic-rich materials

 L = (Mtyp >= 4) & (Mtyp <= 25);

 if any(L)
   Mw(L) = rhoi * phi_i(L) + rhou * phi_u(L);
 end
