 function icXY = regridIC_2D(icA,xA,yA,X,Y,ICmethod)

% Regrids the initial condition values for 2-D simulations from the IC
% input file grid to the model grid.

% Notation:

%   icA         IC values at (xA,yA)        (array)
%   xA          x-locations in input file   (vector)
%   yA          y-locations in input file   (vector)
%   X           X-model grid locations      (vector)
%   Y           Y-model grid locations      (vector)
%   ICmethod    interpolation method        (string)

% Notes:

%   (1) X and Y are generic coordinates.  Thus, this function can be used
%       to regrid icA(rA,zA), icA(xA,yA), or icA(xA,zA).
% ______________________________________________

% Written by:

%   Gary Clow
%   Institute of Arctic and Alpine Research
%   University of Colorado
%   Boulder, Colorado USA
%   Email: gary.clow@colorado.edu
% ______________________________________________

 m  = length(yA);
 N  = length(X);
 M  = length(Y);

 icX  = NaN*ones(m,N);
 icXY = NaN*ones(M,N);

% interpolate along each row

 for k=1:m
   icX(k,:) = interp1(xA,icA(k,:),X,ICmethod);
 end

% interpolate along each column

 for j=1:N
   icXY(:,j) = interp1(yA,icX(:,j),Y,ICmethod);
 end
