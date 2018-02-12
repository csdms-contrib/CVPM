 function rhocp = C_fluid(T,Ftype)

% Finds the volumetric heat capacity of borehole fluids at temperature T.

% Notation:

%   T     = temperature (C)             (vector)
%   Ftype = fluid type                  (scalar)
%   rhocp = volumetric heat capacity    (vector)

% Available fluids:

%   1  water
%   2  JetA
%   3  n-butyl acetate
%   4  Estisol 140
%   5  Estisol 240
%   6  IsoparK
% ______________________________________________

% Written by:

%   Gary Clow
%   Institute of Arctic and Alpine Research
%   University of Colorado
%   Boulder, Colorado USA
%   Email: gary.clow@colorado.edu
% ______________________________________________

 switch Ftype
 case 1                         % water
   rhocp = C_water(T);

 case 2                         % JetA (at -10 C)
   rho   = 829.75;
   cp    = 1817;
   rhocp = rho*cp *ones(size(T));

 case 3                         % n-butyl acetate
   rho   = 901.6 + 0.985 * (0 - T);
   cp    = 2083.6;
   rhocp = rho .* cp;

 case 4                         % Estisol 140
   rho   = 915 + 0.64 * (-30 - T);
   cp    = 2000;
   rhocp = rho*cp *ones(size(T));

 case 5                         % Estisol 240
   rho   = 897 + 0.69 * (-30 - T);
   cp    = 2000;
   rhocp = rho*cp *ones(size(T));

 case 6                         % IsoparK mix
   rho   = 934;
   cp    = 1650;
   rhocp = rho*cp *ones(size(T));
 end
