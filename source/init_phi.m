 function phi = init_phi(C_opt,Z,Mtyp,Marr,col,CS)

% Initializes the porosity.  Except for test materials and ice (Mtyp <=3), 
% the porosity decreases with depth until it reaches a critical value phic.
% Three types of compaction functions are currently available,

%   C_opt = 'off'           no compaction with depth
%         = 'linear'        porosity decreases linearly with depth
%         = 'exponential'   porosity decreases exponentially with depth

% To allow for a variable lithologic profile, different layers are allowed to 
% have different phi0, phic, and hc values.

% Notation:

%   C_opt = compaction function option                  (character string)
%   Z     = depth                                       (M+1)
%   Mtyp  = material type                               (M+1)
%   Marr  = array of material parameters                (M+1,nparams)
%   col   = column within Marr where phi0 is located    (scalar)
%   CS    = coordinate system                           (character string)
%   phi   = porosity                                    (M+1)
%   phi0  = porosity extrapolated to the surface        (M+1)
%   phic  = critical porosity                           (M+1)
%   hc    = compaction length scale                     (M+1)
% ______________________________________________

% Written by:

%   Gary Clow
%   Institute of Arctic and Alpine Research
%   University of Colorado
%   Boulder, Colorado USA
%   Email: gary.clow@colorado.edu
% ______________________________________________

 switch CS
 case 'R'
   Marr = Marr';
 end

% pre-allocate the array assuming zero porosity

 phi = zeros(size(Mtyp));

% reset the porosity for layers consisting of pure ice

 L = (Mtyp >= 2) & (Mtyp <= 3);

 phi(L) = 1;

% reset the porosity for rocks, soils, and organic-rich materials

 L = (Mtyp >= 4) & (Mtyp <= 25);

 switch CS
 case 'R'
   if any(L)
     phi(L) = Marr(L,col);
   end

 case {'Z','RZ','XZ','XYZ'}
   if any(L)

%   unpack phi-related parameters from Marr array

     phi0 = Marr(:,col);   
     phic = Marr(:,col+1);
     hc   = Marr(:,col+2);

% convert to SI base units

     hc = 1000 * hc;        % [m]

%   apply compaction function

     switch C_opt
     case 'off'
       phi(L) = phi0(L);

     case 'linear'
       phi(L) = phi0(L) .* (1 - Z(L) ./ hc(L));

     case 'exponential'
       phi(L) = phi0(L) .* exp(-Z(L) ./ hc(L));
     end

%   apply compaction limit

     L2        = phi < phic;
     L3        = L & L2;
     phi(L3)   = phic(L3);
   end
 end
