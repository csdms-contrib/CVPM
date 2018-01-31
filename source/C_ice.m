 function rhocp = C_ice(T)

% Finds the volumetric heat capacity of ice.

% Notation:

%   T     = temperature (C)                     (vector)
%   rhocp = volumetric heat capacity            (vector)

%   Tr    = reduced temperature (Tk/273.15)
%   Tk    = temperature (K)
% ______________________________________________

% Written by:

%   Gary Clow
%   Institute of Arctic and Alpine Research
%   University of Colorado
%   Boulder, Colorado USA
%   Email: gary.clow@colorado.edu
% ______________________________________________

% parameters

 T0   = 273.15;
 a1   = 2096.1;
 a2   = 1943.8;
 rhoi = 917;

% Tr - 1

 Trm1 = T / T0;

% find the specific heat

 cp = a1 + a2 * Trm1;

% find the volumetric heat capacity

 rhocp    = rhoi * cp;
 L        = T > 0;
 rhocp(L) = 0;
