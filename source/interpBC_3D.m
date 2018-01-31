 function BC = interpBC_3D(t_bcA,bcXY,t,BCmethod)

% Interpolates a boundary condition at location (X,Y) from the tA time grid 
% to time t.  This is for 3-D simulations.
% ______________________________________________

% Written by:

%   Gary Clow
%   Institute of Arctic and Alpine Research
%   University of Colorado
%   Boulder, Colorado USA
%   Email: gary.clow@colorado.edu
% ______________________________________________

% perform interpolation using BCmethod

 BC = interp1(t_bcA,bcXY,t,BCmethod);

% reshape the array

 BC = squeeze(BC);
 BC = BC';
