 function S0 = init_S0(Marr,col,CS)

% Initializes the heat-source parameter S0.

% Notation:

%   Marr = array of material properties                         (M+1,nparams)
%   col  = column within Marr where So is located               (scalar)
%   CS   = coordinate system                                    (character string)
%   S0   = heat-source extrapolated to the surface              (M+1)
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

% extract S0 from the Marr array

 S0 = Marr(:,col);

% convert to SI base units

 S0 = 1e-03 * S0;       % [W/m^3]
