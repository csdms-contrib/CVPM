 function xs0 = init_xs0(Mtyp,Marr,col,CS)

% Initializes the mole fraction of any solutes in the pore spaces if all
% of the pore ice were to melt.

% Notation:

%   Mtyp = material type                            (M+1)
%   Marr = array of material properties             (M+1,nparams)
%   col  = column within Marr where xs0 is located  (scalar)
%   CS   = coordinate system                        (character string)
%   xs0  = mole fraction of solutes with no ice     (M+1)
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

% pre-allocate the array assuming no solutes are present

 xs0 = zeros(size(Mtyp));

% reset the mole fraction for rocks, soils, and organic-rich materials

 L = (Mtyp >= 4) & (Mtyp <= 25);

 if any(L)
   xs0(L) = Marr(L,col);
 end
