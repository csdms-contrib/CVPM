 function Ki = K_ice(T)

% Finds the thermal conductivity of ice at temperature T.
% Based on Yen's (1981) expression.
% Range of validity is 60 K < Tk <= 273.15 K.

% Notation:

%   T  = temperature (C)                (vector)
%   Ki = thermal conductivity of ice    (vector)

% ______________________________________________

% Written by:

%   Gary Clow
%   Institute of Arctic and Alpine Research
%   University of Colorado
%   Boulder, Colorado USA
%   Email: gary.clow@colorado.edu
% ______________________________________________

 Tk = T + 273.15;

% parameters

 a1 = 9.828;
 a2 = 0.0057;

% find the conductivity at temperature T

 Ki    = a1 * exp(-a2 * Tk);
 L     = T > 0;
 Ki(L) = a1 * exp(-a2 * 273.15);

 if any(Tk < 60)
   disp(' ')
   disp('K_ice: temperature is beyond range of validity.')
 end
