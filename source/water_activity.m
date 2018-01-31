 function aw = water_activity(solute,xs)

% Finds the activity of water (aw) for aqueous solutions with various
% electrolytes using the expression provided by Miyawaki et al (1997). 
 
% To prevent the code from crashing, aw values are linearly extrapolated
% when solute concentrations are greater than the max(xs) value for which
% the Miyawaki coefficients are valid.

% Notation:

%   solute  = chemical formula for the solute   (string)
%   xs      = solute mole fraction              (vector)
%   aw      = water activity                    (vector)

% Currently available solutes include,

%   'NaCl'
%   'KCl'
% ______________________________________________

% Written by:

%   Gary Clow
%   Institute of Arctic and Alpine Research
%   University of Colorado
%   Boulder, Colorado USA
%   Email: gary.clow@colorado.edu
% ______________________________________________

 switch solute
 case 'NaCl'
   alpha  =  1.825;
   beta   = -20.78;
   xsmax  =  0.15;      % max(xs) value for alpha,beta validity
   aw_xs  = 0.82565;    % aw @ max(xs)
   dawdxs = -1.6775;    % daw/dxs at max(xs)
 case 'KCl'
   alpha  =  4.754;
   beta   = -49.37;
   xsmax  =  0.07;
   aw_xs  = 0.93593;
   dawdxs = -1.0637;
 otherwise
   disp(' ')
   disp('Water_activity: this solute is currently unavailable.')
   disp(' ')
   pause
 end

% Miyawaki equation

 aw = (1 - xs) .* exp(alpha * xs.^2 + beta * xs.^3);

% use linear extrapolation beyond xsmax

 L     = xs > xsmax;
 aw(L) = aw_xs + dawdxs*(xs(L) - xsmax);

% aw can't be negative

 L     = aw < 0;
 aw(L) = 0;
