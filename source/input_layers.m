 function [zfuL,zfdL,dzL,MtypL,MarrL] = input_layers(Lfile)

% Imports layer information from a comma-delimited *.txt file or from 
% a matlab *.mat file.
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

%   Lfile = Layer file name                                 (character string)
%   zfuL  = upper boundary of layers                        (vector)
%   zfdL  = lower boundary of layers                        (vector)
%   dzL   = vertical resolution requested for each layer    (vector)
%   MtypL = material type                                   (vector)
%   MarrL = material parameters                             (array)
% ______________________________________________

% determine the file type

 ext = findExt(Lfile);

% Extract information from file

 switch ext
 case 'mat'         % matlab mat file

   load(Lfile)

 case 'txt'         % text file

   nc     = 18;                   % number of columns in *.txt file
   fmtstr = repmat('%g,',1,nc);   % format string

   fid = fopen(Lfile,'r');
         fgetl(fid);
         fgetl(fid);
         fgetl(fid);
         fgetl(fid);
   A   = fscanf(fid,fmtstr,[nc inf]);
         fclose(fid);

   A     = A';
   zfuL  = A(:,1);
   zfdL  = A(:,2);
   dzL   = A(:,3);
   MtypL = A(:,4);
   MarrL = A(:,5:nc);
 end
