 function bcXY = regridBC_3D(bcA,xA,yA,X,Y,BCmethod)

% Regrids the boundary condition values for 3-D simulations from the BC 
% input file grid to the model grid.
% ______________________________________________

%	Copyright (C) 2018, Gary Clow

%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, version 3 of the License.

%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License v3.0 for more details.

%   You should have received a copy of the GNU General Public License
%   along with this program.  If not, see <http://www.gnu.org/licenses/>.

%	Developer can be contacted by,
%	email at:
%		gary.clow@colorado.edu
%	paper mail at:
%		Institute of Arctic and Alpine Research
%		University of Colorado
%		Campus Box 450
%		Boulder, CO 80309-0450 USA
% ______________________________________________

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
