 function Kp = K_pore(T,psi_i,psi_u,psi_a,planet)

% Finds the thermal conductivity of the pores within rocks and soils.

% Notation:

%   T     = temperature (C)                             (vector)
%   psi_i = relative volume fraction of ice             (vector)
%   psi_u = relative volume fraction of unfrozen water	(vector)
%   psi_a = relative volume fraction of unfrozen air	(vector)
%   Kp    = conductivity of pores                       (vector)
% ______________________________________________

% Written by:

%   Gary Clow
%   Institute of Arctic and Alpine Research
%   University of Colorado
%   Boulder, Colorado USA
%   Email: gary.clow@colorado.edu
% ______________________________________________

% pre-allocate arrays

 Kp = NaN*ones(size(psi_i));

% find component conductivities

 Ki = K_ice(T);                 % ice
 Ku = K_water(T);               % unfrozen water
 Ka = K_air(T,planet);          % air

% find the conductivity when the pores are filled with ice

 L     = psi_i == 1;
 Kp(L) = Ki(L);

% find the conductivity when the pores are filled with water

 L     = psi_u == 1;
 Kp(L) = Ku(L);

% find the conductivity when the pores are filled with air

 L     = psi_a == 1;
 Kp(L) = Ka(L);

% find the condutivity when there's more than one phase present

 L = (psi_i ~= 1) & (psi_u ~= 1) & (psi_a ~= 1);
 
 Kp(L) = K_pore_3phase(Ki(L),Ku(L),Ka(L),psi_i(L),psi_u(L),psi_a(L));
