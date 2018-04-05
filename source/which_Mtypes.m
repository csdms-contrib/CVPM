 function [Mtypes,m] = which_Mtypes(Mtyp)

% Finds which material types are present.
% ______________________________________________

% Written by:

%   Gary Clow
%   Institute of Arctic and Alpine Research
%   University of Colorado
%   Boulder, Colorado USA
%   Email: gary.clow@colorado.edu
% ______________________________________________

 Mtypes = unique(Mtyp);     % codes for the material types
 m      = length(Mtypes);   % number of unique types present
