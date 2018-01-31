 function bcX = regridBC_2D(bcA,xA,X,BCmethod)

% Regrids the boundary condition values for 2-D simulations from the BC 
% input file grid to the model grid (X,Y,R, or Z).

% Notation:

%   bcA         BC values at (tA,xA)        (m,n)
%   xA          x-locations in input file   (n)
%   X           X-grid locations            (N)
%   BCmethod    interpolation method        (string)
% ______________________________________________

% Written by:

%   Gary Clow
%   Institute of Arctic and Alpine Research
%   University of Colorado
%   Boulder, Colorado USA
%   Email: gary.clow@colorado.edu
% ______________________________________________

 nt = size(bcA,1);
 N  = length(X);

 bcX = NaN*ones(nt,N);

% for each time slice, interpolate bcA onto X-grid

 for i=1:nt
   bcX(i,:) = interp1(xA,bcA(i,:),X,BCmethod);
 end
