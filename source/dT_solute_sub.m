 function theta_s = dT_solute_sub(solute,xs)

% Finds the freezing point depression due to solutes for non-ideal
% solutions using the expression reported by Robinson and Stokes
% (1959).

% Notation:

%   solute  = chemical formula for the solute   (string)
%   xs      = solute mole fraction              (M+1)
%   theta_s = freezing point depression         (M+1)
% ______________________________________________

% Written by:

%   Gary Clow
%   Institute of Arctic and Alpine Research
%   University of Colorado
%   Boulder, Colorado USA
%   Email: gary.clow@colorado.edu
% ______________________________________________

% parameters

 a = 4.76e-06;
 b = 9.687e-3;

% if solutes are present ...

 if ~strcmp(solute,'none')

   switch solute
   case 'NaCl'
     T_eutectic = -21.1;
   otherwise
     disp(' ')
     disp('warning dT_solute: eutectic temperature for this solute is unknown.')
     pause
   end

% find the water activity

   aw = water_activity(solute,xs);
   c  = log(aw);

% find the freezing point depression

   theta_s = (-b + sqrt(b^2 - 4*a*c))/(2*a);

% check for eutectic pt

   L          = theta_s > -T_eutectic;
   theta_s(L) = -T_eutectic;

 else
   theta_s = zeros(size(xs));
 end
