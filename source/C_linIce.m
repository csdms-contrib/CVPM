 function rhocp = C_linIce(T)

% Finds the specific heat and volumetric heat capacity of 'linearized' ice.
% This is the same as C_ice.m except it lacks the phase change at 0 C.

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

 rhocp = rhoi * cp;
