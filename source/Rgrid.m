 function [rf,R,Dr,dr,varep,lambda,Mtyp,Marr,N] ...
    = Rgrid(rfwL,rfeL,drL,MtypL,MarrL,innerBC_type,Rmin,Rmax)

% Sets up the CV grid in the radial dimension.
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

% The layers defined by (rfwL, rfeL) provide boundaries for the control
% volumes.  Each layer can be further subdivided to have a resolution of
% drL meters.

% Layers:

%   rfwL  = inner interface of the layers
%   rfeL  = outer interface of the layers
%   drL   = radial resolution requested for each layer
%   MtypL = type of material
%   MarrL = array of material parameters

% Control volume radial grid:

%   rf     = CV interfaces                      (N+1)
%   R      = CV grid points                     (N+1)
%   Dr     = width of CVs                       (N+1)
%   dr     = distance between grid points       (N+1)
%   lambda = (r_e^2 - r_w^2)/2                  (N+1)
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

% identify layers that are within the (Rmin,Rmax) limits

 jmin = find(rfwL <= Rmin,1,'last');
 if isempty(jmin)
   disp(' ')
   disp('Rgrid error: Rmin < min(rfwL)')
   pause
 end

 jmax = find(rfeL >= Rmax,1,'first');
 if isempty(jmax)
   disp(' ')
   disp('Rgrid error: Rmax > max(rfeL)')
   pause
 end

% set CV interface locations and material parameters at each grid point

 switch innerBC_type
 case {'T','q'}
   rf = [rfwL(jmin); rfwL(jmin)];
 case 'none'
   rf = rfwL(jmin);
 end
 mm = 0;

 for j=jmin:jmax
   Dr  = drL(j);
   rfx = (rfwL(j)+Dr: Dr: min(rfeL(j),Rmax))';
   rf  = [rf; rfx];
   m   = length(rf) - 1;
   Mtyp(mm+1:m,1) = MtypL(j);
   Marr(mm+1:m,:) = repmat(MarrL(j,:),m-mm,1);
   mm  = m;
 end
 clear Dr

 N           = length(rf) -1;
 Mtyp        = Mtyp(1:N+1,1);	% discard unused portion of arrays
 Marr        = Marr(1:N+1,:);
 Mtyp(1,1)   = Mtyp(2,1);       % ensure that values on the boundaries are correctly set
 Marr(1,:)   = Marr(2,:);
 Mtyp(N+1,1) = Mtyp(N,1);
 Marr(N+1,:) = Marr(N,:);

% pre-allocate additional arrays

 R      = NaN*ones(N+1,1);
 Dr     =    zeros(N+1,1);
 dr     =    zeros(N+1,1);
 varep  = NaN*ones(N+1,1);
 lambda = NaN*ones(N+1,1);

% set grid point locations at center of CVs.  Embed a grid point within
% the upper and lower boundaries.

 for k=2:N
   R(k) = 0.5*(rf(k) + rf(k+1));
 end
 R(1)   = rf(1);
 R(N+1) = rf(N+1);

% define width of control volumes

 for k=2:N
   Dr(k) = rf(k+1) - rf(k);
 end

% define distance between grid points

 for k=2:N+1
   dr(k) = R(k) - R(k-1);
 end

% define fractional distance used to compute effective thermal conductivities

 for k=2:N+1
   varep(k) = (R(k) - rf(k)) / dr(k);
 end
 varep(1) = varep(2);

% define radial factor, lambda

 for k=1:N
   lambda(k) = (rf(k+1)^2 - rf(k)^2)/2;
 end

% flip arrays to radial convention

 rf     = rf';
 R      = R';
 Dr     = Dr';
 dr     = dr';
 varep  = varep';
 lambda = lambda';
 Mtyp   = Mtyp';
 Marr   = Marr';
