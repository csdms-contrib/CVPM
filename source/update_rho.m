 function rho = update_rho(Mtyp,rhom,phi,phi_i,phi_u)

% Updates the bulk density field.

% Notation:

%   Mtyp  = material type
%   rhom  = matrix density
%   phi   = porosity
%   phi_i = volume fraction of ice
%   phi_u = volume fraction of unfrozen water

% Notes:
%   (1) We only need to do this for porous materials.
% ______________________________________________

% Written by:

%   Gary Clow
%   Institute of Arctic and Alpine Research
%   University of Colorado
%   Boulder, Colorado USA
%   Email: gary.clow@colorado.edu
% ______________________________________________

% default

 rho = rhom;

% reset bulk density of layers consisting of rocks, soils, or organic-rich layers

 L = (Mtyp >= 4) & (Mtyp <= 25);

 rho(L) = rho_rock(rhom(L),phi(L),phi_i(L),phi_u(L));
