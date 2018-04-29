 function icXY = regridIC_2D(icA,xA,yA,X,Y,ICmethod)

% Regrids the initial condition values for 2-D simulations from the IC
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
