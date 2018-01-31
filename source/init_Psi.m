 function Psi = init_Psi(Mtyp,Marr,col,CS)

% Initializes the relative volume fraction of pores associated with matrix
% particles of radius r.

% Notation:

%   Mtyp = material type                                        (M+1)
%   Marr = array of material properties                         (M+1,nparams)
%   col  = column within Marr where Psi is located              (scalar)
%   CS   = coordinate system                                    (character string)
%   Psi  = volume fraction of pores for particles of radius r   (M+1)
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

% pre-allocate the array assuming no pores are present

 Psi = zeros(size(Mtyp));

% reset the mole fraction for rocks, soils, and organic-rich materials

 L = (Mtyp >= 4) & (Mtyp <= 25);

 if any(L)
   Psi(L) = Marr(L,col);
 end
