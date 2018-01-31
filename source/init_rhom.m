 function rhom = init_rhom(Mtyp,Marr,col,CS)

% Initializes the density of the matrix material.

% Notation:

%   Mtyp = material type                            (M+1)
%   Marr = array of material parameters             (M+1,nparams)
%   col  = column within Marr where rhm0 is located (scalar)
%   CS   = coordinate system                        (character string)
%   rhom = matrix density                           (M+1)

% Notes:
%   (1) For ice (Mtyp=3), ignore what's found in the *_layers.txt file.
% ______________________________________________

% Written by:

%   Gary Clow
%   Institute of Arctic and Alpine Research
%   University of Colorado
%   Boulder, Colorado USA
%   Email: gary.clow@colorado.edu
% ______________________________________________

% set rhom according to what's found in the *_layers.txt file

 switch CS
 case 'R'
   rhom = Marr(col,:);
 otherwise
   rhom = Marr(:,col);
 end

% reset rhom if layer consists of pure ice

 L = (Mtyp >= 2) & (Mtyp <= 3);

 if any(L)
   rhom(L) = 917;
 end
