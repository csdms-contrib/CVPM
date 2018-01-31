 function theta_s = dT_solute(Mtyp,solute,xs0,phi_w,phi_u)

% Finds the freezing point depression due to solutes.

% Notation:

%   Mtyp    = material type                             (M+1)
%   solute  = chemical formula for the solute           (string)
%   xs0     = solute mole fraction with no ice          (M+1)
%   phi_w   = volume fraction of water (solid+liquid)   (M+1)
%   phi_u   = volume fraction of unfrozen water         (M+1)
%   theta_s = freezing point depression                 (M+1)

 dflag = 0;     % (0) quiet mode, (1) dumps various fields to screen.
% ______________________________________________

% Written by:

%   Gary Clow
%   Institute of Arctic and Alpine Research
%   University of Colorado
%   Boulder, Colorado USA
%   Email: gary.clow@colorado.edu
% ______________________________________________

% pre-allocate arrays

 theta_s = zeros(size(Mtyp));

% find the solute mole fraction

 xs = xs0 .* (phi_w ./ phi_u);

% find the freezing point depression

 L = (Mtyp >=4) & (Mtyp <= 25);

 if any(L)
   theta_s(L) = dT_solute_sub(solute,xs(L));
 end

 if dflag
   disp(' ')
   disp('dT_solute: Mtyp phi_w phi_u xs theta_s')
   [Mtyp phi_w phi_u xs theta_s]
   pause(2)
 end
