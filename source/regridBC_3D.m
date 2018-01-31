 function bcXY = regridBC_3D(bcA,xA,yA,X,Y,BCmethod)

% Regrids the boundary condition values for 3-D simulations from the BC 
% input file grid to the model grid.

% Notation:

%   bcA         BC values at (tA,xA,yA)     (array)
%   xA          x-locations in input file   (vector)
%   yA          y-locations in input file   (vector)
%   X           X-model grid locations      (vector)
%   Y           Y-model grid locations      (vector)
%   BCmethod    interpolation method        (string)

% Notes:

%   (1) X and Y are generic coordinates.  Thus, this function can be used
%       to regrid bcA(xA,yA), bcA(xA,zA), or bcA(yA,zA).
% ______________________________________________

% Written by:

%   Gary Clow
%   Institute of Arctic and Alpine Research
%   University of Colorado
%   Boulder, Colorado USA
%   Email: gary.clow@colorado.edu
% ______________________________________________

 nt = size(bcA,1);
 m  = length(yA);
 N  = length(X);
 M  = length(Y);

 bcX  = NaN*ones(nt,m,N);
 bcXY = NaN*ones(nt,M,N);

% For each time slice, interpolate bcA onto X-grid

 for i=1:nt

% interpolate along each row

   for k=1:m
     bcX(i,k,:) = interp1(xA,squeeze(bcA(i,k,:)),X,BCmethod);
   end

% interpolate along each column

   for j=1:N
     bcXY(i,:,j) = interp1(yA,squeeze(bcX(i,:,j)),Y,BCmethod);
   end
 end
