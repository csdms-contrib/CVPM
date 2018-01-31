 function K = K_rock_porous(T,Km,phi,phi_i,phi_u,planet)

% Finds the bulk thermal conductivity of porous rocks and soils at 
% temperature T.

% Notation:

%   T     = temperature (C)                         (vector)
%   Km    = mineral grain conductivity              (vector)
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

% find the relative volume fractions

 psi_i    = phi_i ./ phi;           % ice
 psi_u    = phi_u ./ phi;           % unfrozen water
 psi_a    = 1 - (psi_i + psi_u);    % air
 L        = psi_a < 0;
 psi_a(L) = 0;                      % make sure psi_a doesn't go slightly negative due to numerics

% find the conductivity of the pores

 Kp = K_pore(T,psi_i,psi_u,psi_a,planet);

% find the bulk conductivity of the mineral-matrix/pore mixture

 K = K_BM2r(Km,Kp,phi);
