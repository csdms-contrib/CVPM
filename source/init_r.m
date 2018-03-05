 function r = init_r(Mtyp,Marr,col,CS)

% Initializes the effective radius of the matrix particles.

% Notation:

%   Mtyp = material type                            (M+1)
%   Marr = array of material properties             (M+1,nparams)
%   col  = column within Marr where r is located	(scalar)
%   CS   = coordinate system                        (character string)
%   r    = particle radius                          (M+1)
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

% pre-allocate diameters assuming no pores are present

 d = zeros(size(Mtyp));

% reset the radius for rocks, soils, and organic-rich materials

 L = (Mtyp >= 4) & (Mtyp <= 25);

 if any(L)
   d(L) = Marr(L,col);
 end

 r = d/2;

% convert to SI base units

 r = 1e-06 * r;     % [m]
