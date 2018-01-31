 function [rhocp,C] = C_veg(T,Mtyp,rhom,cpm0,phi,phi_i,phi_u,dphiudT)

% Finds the volumetric heat capacity of organic-rich material at 
% temperature T.

% Notation:

%   T       = temperature (C)                           (vector)
%   Mtyp    = material type                             (vector)
%   rhom    = mean matrix density                       (vector)
%   cpm0    = pecific heat of matrix at 20 C            (vector)
%   phi     = porosity                                  (vector)
%   phi_i   = volume fraction of ice                    (vector)
%   phi_u   = volume fraction of unfrozen water         (vector)
%   dphiudT = d(phi_u)/dT                               (vector)
%   rhocp   = lattice component of the heat capacity	(vector)
%   C       = volumetric heat capacity                  (vector)
% ______________________________________________

% Written by:

%   Gary Clow
%   Institute of Arctic and Alpine Research
%   University of Colorado
%   Boulder, Colorado USA
%   Email: gary.clow@colorado.edu
% ______________________________________________

% parameters

 rhow = 1000;           % density of water
 Lf   = 3.34e05;        % latent heat of fusion (water)

% find component heat capacities

 rhocpm = C_organic(T,Mtyp,rhom,cpm0);  % organic-rich material
 rhocpi = C_ice(T);                     % ice
 rhocpu = C_water(T);                   % unfrozen water

% find the lattice-vibration component of the heat capacity

 rhocp = (1 - phi).*rhocpm + (phi_i).*rhocpi + (phi_u).*rhocpu;

% find the phase change contribution

 Cphase = rhow * Lf * dphiudT;

% find the total volumetric heat capacity

 C = rhocp  + Cphase;
