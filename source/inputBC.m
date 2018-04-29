 function [descript,t,t_units,BC,Imethod,Sc1,Sc2] = inputBC(BCfile,CS)

% Imports a time-dependent boundary condition.  
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

% Acceptable file types are:

%   (1) a comma delimited *.txt file (1-D, 2-D simulations)
%   (2) a matlab *.mat file

% Notation:

%   BCfile   = BC file name                                 (character string)
%   CS       = coordinate system                            (character string)

%   descript = boundary condition descriptor                (character string)
%   t        = time                                         (vector)
%   t_units  = time units                                   (character string)
%   BC       = BC values                                    (array)
%   Imethod  = method to be used when interpolating the BC  (character string)
%   Sc1      = values of 1st space coordinate returned      (vector)
%   Sc2      = values of 2nd space coordinate returned      (vector)

% Notes:

%   (1) For 1-D simulations, Sc1 and Sc2 are null.
%   (2) For 2-D simulations, the BC is specified as a function of R,X, or Z.
%       This coordinate is loaded into Sc1 while Sc2 is null.
%   (3) For 3-D simulations, the BC is specified as a function XZ,YZ, or XY.
%       The respective coordinate vectors are loaded into Sc1 and Sc2.
% ______________________________________________

% determine the file type

 ext = findExt(BCfile);

% Extract information from file

 switch ext
 case 'mat'         % *.mat file ------------------------

   load(BCfile)

% assign space coordinates
   switch CS
   case {'R','X','Z'}
     Sc1 = [];
     Sc2 = [];
   case {'RZ','XZ'}
     if exist('R','var')
       Sc1 = R;
     elseif exist('X','var')
       Sc1 = X;
     elseif exist('Y','var')
       Sc1 = Y;
     elseif exist('Z','var')
       Sc1 = Z;
     end
     Sc2 = [];
   case 'XYZ'
     if exist('X','var')
       Sc1 = X;
       if exist('Y','var')
         Sc2 = Y;
       else
         Sc2 = Z;
       end
     else
       Sc1 = Y;
       Sc2 = Z;
     end
   end

 case 'txt'         % text file --------------------------

   fid = fopen(BCfile,'r');

% read descriptors

   descript1 = fgetl(fid);
   descript2 = fgetl(fid);
               fgetl(fid);
   descript  = [descript1 '; ' descript2];

% time units

   s = fgetl(fid); S = extract_strs(s); t_units = S{1};

% interpolation method

   s = fgetl(fid); S = extract_strs(s); Imethod = S{1};

% assign space coordinates

   switch CS
   case {'R','X','Z'}
     Sc1 = [];
     Sc2 = [];
   case {'RZ','XZ'}
     s = fgetl(fid); Sc1 = extract_nums(s);
     Sc2 = [];
   case 'XYZ'
     disp(' ')
     disp('input_BC: for 3D coordinate systems, use make3D_BC to create the BC file.')
     pause
   end

% read in BC, one time record (line) at a time

   A = [];
   while feof(fid) == 0
     s = fgetl(fid);
     x = extract_nums(s);
     A = [A; x];
   end
   n  = size(A,2);
   t  = A(:,1);
   BC = A(:,2:n);
 end
