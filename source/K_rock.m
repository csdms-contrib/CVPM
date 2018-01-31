 function K = K_rock(T,Mtyp,Km0,phi,phi_i,phi_u,planet)

% Finds the bulk thermal conductivity of rocks and soils at temperature T.

% Notation:

%   T     = temperature (C)                         (vector)
%   Mtyp  = material type                           (vector)
%   Km0   = matrix conductivity at 0 C              (vector)
%   phi   = porosity                                (vector)
%   phi_i = volume fraction of ice                  (vector)
%   phi_u = volume fraction of unfrozen water       (vector)
%   K     = bulk thermal conductivity               (vector)
% ______________________________________________

% Written by:

%   Gary Clow
%   Institute of Arctic and Alpine Research
%   University of Colorado
%   Boulder, Colorado USA
%   Email: gary.clow@colorado.edu
% ______________________________________________

% pre-allocate arrays

 K = NaN*ones(size(phi));

% find the conductivity of the mineral grains

 Km = K_mineral(T,Mtyp,Km0);

% find the bulk conductivity of non-porous rocks

 L    = phi == 0;
 K(L) = Km(L);

% find the bulk conductivity of porous rocks and soils

 L    = phi ~= 0;
 K(L) = K_rock_porous(T(L),Km(L),phi(L),phi_i(L),phi_u(L),planet);
