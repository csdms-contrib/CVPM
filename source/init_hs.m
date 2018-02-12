 function hs = init_hs(Marr,col,CS)

% Initializes the heat-source length scale.

% Notation:

%   Marr = array of material properties             (M+1,nparams)
%   col  = column within Marr where So is located	(scalar)
%   CS   = coordinate system                        (character string)
%   hs   = heat-source length scale              	(M+1)
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

% extract hs from the Marr array

 hs = Marr(:,col);

% convert to SI base units

 hs = 1000 * hs;        % [m]
