 function [xf,X,Dx,dx,varep,Mtyp,Marr,N] ...
    = Xgrid(xfwL,xfeL,dxL,MtypL,MarrL,Xmin,Xmax)

% Sets up the CV grid in a horizontal dimension.  It can be used to setup 
% either the X- or Y-grid for the cartesian coordinate system.
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

% The layers defined by (xfwL, xfeL) provide boundaries for the control
% volumes.  Each layer can be further subdivided to have a resolution of
% dxL meters.

% Layers:

%   xfwL  = inner interface of the layers
%   xfeL  = outer interface of the layers
%   dxL   = horizontal resolution requested for each layer
%   MtypL = type of material
%   MarrL = array of material parameters

% Control volume horizontal x-grid:

%   xf     = CV interfaces                      (N+1)
%   X      = CV grid points                     (N+1)
%   Dx     = width of CVs                       (N+1)
%   dx     = distance between grid points       (N+1)
%   varep  = fractional distance                (N+1)
%   Mtyp   = material type at each grid point	(N+1)
%   Marr   = array of material parameters       (nparams,N+1)
%   nparams = number of material properties     (scalar)

% Notes:

%   (1) There are (N-1) CVs and (N+1) CV interfaces.
%   (2) Full control volumes are used adjacent to both boundaries.
%   (3) A grid point is embedded within both the inner and outer boundaries.
%   (4) Arrays keyed to the control volumes (Mtyp,Marr) have values from
%       the adjacent CVs stored in their first and last elements, eg.
%       Mtyp(1) = Mtyp(2) and Mtyp(N+1) = Mtyp(N).
% ______________________________________________

% pre-allocate arrays

 nparams = size(MarrL,2);       % number of material parameters
 Mtyp    = NaN*ones(1000,1);
 Marr    = NaN*ones(1000,nparams);

% identify layers that are within the (Xmin,Xmax) limits

 jmin = find(xfwL <= Xmin,1,'last');
 if isempty(jmin)
   disp(' ')
   disp('Xgrid error: Xmin < min(xfwL)')
   pause
 end

 jmax = find(xfeL >= Xmax,1,'first');
 if isempty(jmax)
   disp(' ')
   disp('Xgrid error: Xmax > max(xfeL)')
   pause
 end

% set CV interface locations and material parameters at each grid point

 xf = [xfwL(1); xfwL(1)];
 mm = 0;

 for j=jmin:jmax
   Dx  = dxL(j);
   xfx = (xfwL(j)+Dx: Dx: min(xfeL(j),Xmax))';
   xf  = [xf; xfx];
   m   = length(xf) - 1;
   Mtyp(mm+1:m,1) = MtypL(j);
   Marr(mm+1:m,:) = repmat(MarrL(j,:),m-mm,1);
   mm  = m;
 end
 clear Dx

 N           = length(xf) -1;
 Mtyp        = Mtyp(1:N+1,1);	% discard unused portion of arrays
 Marr        = Marr(1:N+1,:);
 Mtyp(1,1)   = Mtyp(2,1);       % ensure that values on the boundaries are correctly set
 Marr(1,:)   = Marr(2,:);
 Mtyp(N+1,1) = Mtyp(N,1);
 Marr(N+1,:) = Marr(N,:);

% pre-allocate additional arrays

 X      = NaN*ones(N+1,1);
 Dx     =    zeros(N+1,1);
 dx     =    zeros(N+1,1);
 varep  = NaN*ones(N+1,1);

% set grid point locations at center of CVs.  Embed a grid point within
% the upper and lower boundaries.

 for k=2:N
   X(k) = 0.5*(xf(k) + xf(k+1));
 end
 X(1)   = xf(1);
 X(N+1) = xf(N+1);

% define width of control volumes

 for k=2:N
   Dx(k) = xf(k+1) - xf(k);
 end

% define distance between grid points

 for k=2:N+1
   dx(k) = X(k) - X(k-1);
 end

% define fractional distance used to compute effective thermal conductivities

 for k=2:N+1
   varep(k) = (X(k) - xf(k)) / dx(k);
 end
 varep(1) = varep(2);

% flip arrays to horizontal convention

 xf     = xf';
 X      = X';
 Dx     = Dx';
 dx     = dx';
 varep  = varep';
 Mtyp   = Mtyp';
 Marr   = Marr';
