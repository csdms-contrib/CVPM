 function Sr = init_Sr(Mtyp,Marr,col,CS)

% Initializes the initial degree of water (both phases) saturation.

% Notation:

%   Mtyp = material type                            (M+1)
%   Marr = array of material properties             (M+1,nparams)
%   col  = column within Marr where Sr is located   (scalar)
%   CS   = coordinate system                        (character string)
%   Sr   = degree of saturation                     (M+1)
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

% pre-allocate the array assuming complete saturation

 Sr = ones(size(Mtyp));

% reset the saturation for rocks, soils, and organic-rich materials

 L = (Mtyp >=4) & (Mtyp <= 25);

 if any(L)
   Sr(L) = Marr(L,col);
 end
