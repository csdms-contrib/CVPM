 function [zf,Z,Dz,dz,varep,Mtyp,Marr,M] ...
      = Zgrid(zfuL,zfdL,dzL,MtypL,MarrL,Zmin,Zmax)

% Sets up the CV grid in the vertical dimension.

% The layers defined by (zfuL, zfdL) provide boundaries for the control
% volumes.  Each layer can be further subdivided to have a resolution of
% dzL meters.

% Layers:

%   zfuL  = upper interface of the layers
%   zfdL  = lower interface of the layers
%   dzL   = vertical resolution requested for each layer
%   MtypL = type of material
%   MarrL = array of material parameters

% Control volume vertical grid:

%   zf    = CV interfaces                       (M+1)
%   Z     = CV grid points                      (M+1)
%   Dz    = width of CVs                        (M+1)
%   dz    = distance between grid points        (M+1)
%   varep = fractional distance                 (M+1)
%   Mtyp  = material type at each grid point    (M+1)
%   Marr  = array of material parameters        (M+1,nparams)
%   nparams = number of material parameters     (scalar)

% Notes: 

%   (1) There are (M-1) CVs and (M+1) CV interfaces.
%   (2) Full control volumes are used adjacent to both boundaries.
%   (3) A grid point is embedded within both the upper and lower boundaries.
%   (4) Arrays keyed to the control volumes (Mtyp,Marr) have values from
%       the adjacent CVs stored in their first and last elements, eg.
%       Mtyp(1) = Mtyp(2) and Mtyp(M+1) = Mtyp(M).
% ______________________________________________

% Written by:

%   Gary Clow
%   Institute of Arctic and Alpine Research
%   University of Colorado
%   Boulder, Colorado USA
%   Email: gary.clow@colorado.edu
% ______________________________________________

% pre-allocate arrays

 nparams = size(MarrL,2);
 Mtyp    = NaN*ones(1000,1);
 Marr    = NaN*ones(1000,nparams);

% identify layers that are within the (Zmin,Zmax) limits

 jmin = find(zfuL <= Zmin,1,'last');
 if isempty(jmin)
   disp(' ')
   disp('Zgrid error: Zmin < min(zfuL)')
   pause
 end

 jmax = find(zfdL >= Zmax,1,'first');
 if isempty(jmax)
   disp(' ')
   disp('Zgrid error: Zmax > max(zfdL)')
   pause
 end

% set CV interface locations and material parameters at each grid point

 zf = [zfuL(1); zfuL(1)];
 mm = 0;

 for j=jmin:jmax
   Dz  = dzL(j);
   zfx = (zfuL(j)+Dz: Dz: min(zfdL(j),Zmax))';
   zf  = [zf; zfx];
   m   = length(zf) - 1;
   Mtyp(mm+1:m,1) = MtypL(j);
   Marr(mm+1:m,:) = repmat(MarrL(j,:),m-mm,1);
   mm  = m;
 end
 clear Dz

 M           = length(zf) -1;
 Mtyp        = Mtyp(1:M+1,1);	% discard unused portion of arrays
 Marr        = Marr(1:M+1,:);
 Mtyp(1,1)   = Mtyp(2,1);       % ensure that values on the boundaries are correctly set
 Marr(1,:)   = Marr(2,:);
 Mtyp(M+1,1) = Mtyp(M,1);
 Marr(M+1,:) = Marr(M,:);

% pre-allocate additional arrays

 Z     = NaN*ones(M+1,1);
 Dz    =    zeros(M+1,1);
 dz    =    zeros(M+1,1);
 varep = NaN*ones(M+1,1);

% set grid point locations at center of CVs.  Embed a grid point within
% the upper and lower boundaries.

 for k=2:M
   Z(k) = 0.5*(zf(k) + zf(k+1));
 end
 Z(1)   = zf(1);
 Z(M+1) = zf(M+1);

% define width of control volumes

 for k=2:M
   Dz(k) = zf(k+1) - zf(k);
 end

% define distance between grid points

 for k=2:M+1
   dz(k) = Z(k) - Z(k-1);
 end

% define fractional distance used to compute effective thermal conductivities

 for k=2:M+1
   varep(k) = (Z(k) - zf(k)) / dz(k);
 end
 varep(1) = varep(2);
