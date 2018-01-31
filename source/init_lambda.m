 function lambda = init_lambda(Mtyp,Marr,col,CS)

% Initializes the interfacial melting parameter.

% Notation:

%   Mtyp = material type                                (M+1)
%   Marr = array of material properties                 (M+1,nparams)
%   col  = column within Marr where lambda is located   (scalar)
%   CS   = coordinate system                            (character string)
%   lambda = interfacial melting parameter              (M+1)
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

 lambda = zeros(size(Mtyp));

% reset the interfacial melting parameter for rocks, soils, and organic-rich materials

 L = (Mtyp >= 4) & (Mtyp <= 25);

 if any(L)
   lambda(L) = Marr(L,col);
 end

% convert to SI base units

 lambda = 1e-06 * lambda;       % [m K^{1/3}]
