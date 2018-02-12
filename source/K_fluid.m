 function K = K_fluid(T,Ftype)

% Finds the thermal conductivity of borehole fluids at temperature T.

% Notation:

%   T     = temperature (C)          (vector)
%   Ftype = fluid type               (scalar)
%   K     = thermal conductivity     (vector)

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
 case 1                     % water
   K = K_water(T);

 case 2                     % JetA (at -10 C)
   K = 0.12046*ones(size(T));

 case 3                     % n-butyl acetate (at -10 C)
   K  = 0.1434*ones(size(T));

 case {4,5}                 % Estisol 140/240
   K = 0.18*ones(size(T));

 case 6                     % IsoparK mix
   K = 0.15*ones(size(T));
 end
