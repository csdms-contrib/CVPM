 function cpm0 = init_cpm(Mtyp,Marr,col,CS)

% Initializes the specific heat of the matrix material at reference 
% temperature (20 C).

% Notation:

%   Mtyp = material type                            (M+1)
%   Marr = array of material parameters             (M+1,nparams)
%   col  = column within Marr where cpm0 is located (scalar)
%   CS   = coordinate system                        (character string)
%   cpm0 = matrix specific heat at 20 C             (M+1)

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

% set cpm0 according to what's found in the *_layers.txt file

 switch CS
 case 'R'
   cpm0 = Marr(col,:);
 otherwise
   cpm0 = Marr(:,col);
 end

% reset cmp0 if layer consists of pure ice

 L = Mtyp == 3;

 if any(L)
   cpm0(L) = NaN;
 end
