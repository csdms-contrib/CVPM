 function rho = init_rho(Mtyp,rhom,phi,phi_i,phi_u)

% Initializes the bulk density.

% Notation:

%   Mtyp  = material type                               (M+1)
%   rhom  = density of matrix material                  (M+1)
%   phi   = porosity                                    (M+1)
%   phi_i = volume fraction of ice                      (M+1)
%   phi_u = volume fraction of unfrozen water           (M+1)
%   rho   = bulk density                                (M+1)

% Notes:

%   (1) This is an estimate of rho based on the assumption that all the
%       water is in the ice phase.  When init_rho.m is called, we really
%       don't know how much water is ice and how much is liquid.
% ______________________________________________

% Written by:

%   Gary Clow
%   Institute of Arctic and Alpine Research
%   University of Colorado
%   Boulder, Colorado USA
%   Email: gary.clow@colorado.edu
% ______________________________________________

% parameters

 rhoi =  917;       % density of ice
 rhou = 1000;       % unfrozen water

% pre-allocate the array to the matrix density

 rho = rhom;

% reset bulk density for rocks, soils, and organic-rich layers

 L = (Mtyp >= 4) & (Mtyp <= 25);

 if any(L)

   rho(L) = (1 - phi(L)).* rhom(L) + phi_i(L)*rhoi + phi_u(L)*rhou;
 end
