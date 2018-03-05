 function n21 = init_n21(Mtyp,Marr,col,CS)

% Initializes (n2/n1), the ratio of the number of particles (or pores) with
% radius r2 to those with radius r1.

% Notation:

%   Mtyp = material type                                (M+1)
%   Marr = array of material properties                 (M+1,nparams)
%   col  = column within Marr where n21 is located      (scalar)
%   CS   = coordinate system                            (character string)
%   n21  = (n2/n1)                                      (M+1)
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

 n21 = zeros(size(Mtyp));

% reset n21 for rocks, soils, and organic-rich materials

 L = (Mtyp >= 4) & (Mtyp <= 25);

 if any(L)
   n21(L) = Marr(L,col);
 end
