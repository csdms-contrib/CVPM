% CVPM_view.m

% View the results from a CVPM simulation.
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

% pickup the location of the working directory from the CVPM.config file

 [wdir,~,~,~] = input_config;
 addpath([wdir '/source'],'-end')
 addpath([wdir '/utilities'],'-end')
 cd(wdir)
% ______________________________________________

% Select coordinate system

 disp(' ')
 disp('CVPM coordinate systems: ')
 disp('[1] 1D radial      (R)')
 disp('[2] 1D vertical    (Z)')
 disp('[3] 2D cylindrical (RZ)')
 disp('[4] 2D cartesian   (XZ)')
 disp('[5] 3D catesian    (XYZ)')
 CSopt = input('Select CS option number (1,2,...): ');

% Launch the appropriate viewer

 switch CSopt
 case 1
   view_R
 case 2
   view_Z
 case 3
   view_RZ
 case 4
   view_XZ
 case 5
   view_XYZ
 end
