 function Ki = K_linIce(T,K0)

% Finds the thermal conductivity of a substance for which K is linearly
% dependent on temperature.  For this we choose properties similar to ice.

% Notation:

%   T  = temperature (C)                      (vector)
%   K0 = conductivity @ 0 C                   (vector)
%   Ki = conductivity at temperature T        (vector)
% ______________________________________________

% Written by:

%   Gary Clow
%   Institute of Arctic and Alpine Research
%   University of Colorado
%   Boulder, Colorado USA
%   Email: gary.clow@colorado.edu
% ______________________________________________

% parameters

 T0   = 0;              % reference temperature (C)
 dKdT = -0.0125;        % dK/dT

% find the conductivity

 Ki = K0 + dKdT *(T - T0);
