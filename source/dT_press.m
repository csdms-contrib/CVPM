 function theta_p = dT_press(planet,P_opt,Z,dz,rho)

% Finds the freezing point depression due to pressure.  Pore pressures
% can be calculated assuming either 'hydrostatic' or 'lithostatic' 
% environmental conditions.

% Notation:

%   planet                                   (string)
%   P_opt   = pressure option                (string)
%   Z       = depth                          (M+1)
%   dz      = thickness of each CV           (M+1)
%   rho     = bulk density                   (M+1)
%   theta_p = feezing point depression       (M+1)

% Notes:

% (1) planet = 'earth' or 'mars'
% (2) P_opt  = 'hydrostatic' or 'lithostatic'
% ______________________________________________

% Written by:

%   Gary Clow
%   Institute of Arctic and Alpine Research
%   University of Colorado
%   Boulder, Colorado USA
%   Email: gary.clow@colorado.edu
% ______________________________________________

% parameters

 a    = 9.8e-08;        % air-saturated pressure sensitivity (K/Pa)
 P0   = 611.657;        % triple point pressure (Pa)
 rhow = 1000;           % density of water (kg/m^3)

 switch planet
 case 'earth'
   g    = 9.78;         % surface gravity (m/s^2)
   Patm = 1.01e05;      % surface air pressure (Pa)
 case 'mars'
   g    = 3.72;
   Patm = 700;
 otherwise
   disp(' ')
   disp('warning dT_press: unknown planet.')
   pause
 end

% find pore pressure

 switch P_opt
 case 'off'
   P = P0 *ones(size(Z));
 case 'hydrostatic'
   rhoZ = rhow * Z;
   P    = Patm + g * rhoZ;
 case 'lithostatic'
   rhodz = rho .* dz;
   rhoZ  = zeros(size(Z));
   for i=1:length(Z)
     rhoZ(i) = sum(rhodz(1:i));
   end
   P = Patm + g * rhoZ;
 otherwise
   disp(' ')
   disp('warning dT_press: unknown pressure option.')
   pause
 end

% find freezing point depression

 theta_p = a * (P - P0);
