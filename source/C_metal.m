 function rhocp = C_metal(T,Mtyp)

% Finds the volumetric heat capacity of metals.

% Notation:

%   T     = temperature (C)             (vector)
%   Mtyp  = meterial type               (scalar)
%   rhocp = volumetric heat capacity    (vector)

% Available metals:

%   40  steel drill pipe
%   41  stainless steel
%   42  cast iron
%   43  aluminum
%   44  copper

% Notes:
%   (1) This routine does not currently take into account the temperature 
%       dependence.

% Source: internet tables
% ______________________________________________

% Written by:

%   Gary Clow
%   Institute of Arctic and Alpine Research
%   University of Colorado
%   Boulder, Colorado USA
%   Email: gary.clow@colorado.edu
% ______________________________________________

 switch Mtyp(1)
 case 40            % steel drill pipe
   rho = 7833;
   cp  = 490;

 case 41            % stainless steel
   rho = 7982;
   cp  = 480;

 case 42            % cast iron
   rho = 7300;
   cp  = 460;

 case 43            % aluminum
   rho = 2705;
   cp  = 910;

 case 44            % copper
   rho = 8944;
   cp  = 390;

 end
 rhocp = rho*cp *ones(size(T));
