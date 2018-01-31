 function rho = rho_rock(rhom,phi,phi_i,phi_u)

% Finds the density of porous rocks and soils.

% Notation:

%   rhom  = matrix density                          (vector)
%   phi   = porosity                                (vector)
%   phi_i = volume fraction of ice                  (vector)
%   phi_u = volume fraction of unfrozen water       (vector)
%   rho   = bulk density                            (vector)
% ______________________________________________

% Written by:

%   Gary Clow
%   Institute of Arctic and Alpine Research
%   University of Colorado
%   Boulder, Colorado USA
%   Email: gary.clow@colorado.edu
% ______________________________________________

% component densities

 rhoi =  917;       % ice
 rhou = 1000;       % unfrozen water

% find bulk density

 rho = (1 - phi).*rhom + (phi_i)*rhoi + (phi_u)*rhou;
