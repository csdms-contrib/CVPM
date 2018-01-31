 function K = K_metal(T,Mtyp)

% Finds the thermal conductivity of metals.

% Notation:

%   T    = temperature (C)          (vector)
%   Mtyp = metal type               (scalar)
%   K    = thermal conductivity     (vector)

% Available metals:

%   40  steel drill pipe
%   41  stainless steel
%   42  cast iron
%   43  aluminum
%   44  copper

% Notes:
%   (1) This routine does not currently take into account the temperature 
%       dependence.

% Source: wikipedia (for all values except steel drill pipe)
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
   K = 45;

 case 41            % stainless steel (note: values vary considerably)
   K = 16.5;

 case 42            % cast iron
   K = 55;

 case 43            % aluminum
   K = 236;

 case 44            % copper
   K = 401;
 end

 K = K * ones(size(Mtyp));
