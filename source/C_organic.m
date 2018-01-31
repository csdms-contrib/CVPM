 function rhocp = C_organic(T,Mtyp,rhom,cpm0)

% Finds the volumetric heat capacity of a organic-rich matrix at temperature
% T.  The matrix is assumed to consist of a mixture of peat and mineral
% particles.

% Notation:

%   T     = temperature (C)                     (vector)
%   Mtyp  = material type                       (vector)
%   rhom  = mean density of mineral grains      (vector)
%   cpm0  = specific heat of minerals (20 C)    (vector)
%   rhocp = volumetric heat capacity            (vector)

% Available matrix types (Mtyp):

%   20      100% peat
%   21  	 75% peat/25% mineral mix
%   22       50% peat/50% mineral mix
%   23       25% peat/75% mineral mix
% ______________________________________________

% Written by:

%   Gary Clow
%   Institute of Arctic and Alpine Research
%   University of Colorado
%   Boulder, Colorado USA
%   Email: gary.clow@colorado.edu
% ______________________________________________

 Tk = T + 273.15;       % kelvin temperature

% heat capacity of peat particles (see 2004Ling_a).

 rhocpo = 1e06*(-0.1333 + 6.255e-03 * Tk);

% heat capacity of mineral particles

 rhocpm = C_mineral(T,Mtyp,rhom,cpm0);

% heat capacity of mixture

 switch Mtyp(1)
 case 20
   f = 1;
 case 21
   f = 0.75;
 case 22
   f = 0.5;
 case 23
   f = 0.25;
 case 24
   f = 0;
 end

 rhocp = f *rhocpo + (1-f) *rhocpm;
