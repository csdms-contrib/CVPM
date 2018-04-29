 function [wdir,CS,opt,experimC] = input_config

% Imports information from the CVPM config file.
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

%   wdir    = working directory             (character string)
%   CS      = coordinate system             (character string)
%   opt     = configuration options         (vector)
%   experim = name of simulation            (character string)
% ______________________________________________

 Cfile = 'CVPM.config';
 fid   = fopen(Cfile,'r');

% read past headers
 fgetl(fid);
 fgetl(fid);

% working directory
 s = fgetl(fid); S = extract_strs(s); wdir = S{1};

% coordinate_system
 s = fgetl(fid); S = extract_strs(s); CS = S{1};

% configuration options
 s = fgetl(fid); opt = extract_nums(s);

% experiment(s)

 ic = 1;    % experiment counter

 while feof(fid) == 0
   s  = fgetl(fid); S = extract_strs(s); experimC{ic} = S{1};
   ic = ic + 1;
 end

 fclose(fid);
