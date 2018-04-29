% CVPM_portal.m

% This is the portal to the CVPM permafrost modeling system.
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

% This script will run CPS, CVPM, or both with a single call.
% Multiple experiments can be run depending on what's found in the 
% CVPM.config file.

% Options:

% Gopt = 1  Zlayers file is in 'wdir/geo/', text format
%        2  Zlayers file is in 'wdir/geo/', matlab format
%        3  Zlayers file is in 'wdir/tmp/', matlab format

% Ropt = 1  run CPS on one or more files
%        2  run CVPM on one or more files
%        3  run CPS and CVPM on one or more files
% ______________________________________________

% set working directory and add paths to required subdirectories

 [wdir,CS,opt,experimC] = input_config;
 addpath( wdir,                '-end')
 addpath([wdir '/source'],'-end')
 addpath([wdir '/utilities'],'-end')
 savepath
 cd(wdir)

 disp(' ')
 disp(['Working Directory: ' wdir])
 disp(' ')
 disp(['Coordinate System: ' CS])

% set run options

 Gopt = opt(1);
 Ropt = opt(2);

% launch CPS/CVPM for each experiment

 nexperims = size(experimC,2);

 for i=1:nexperims
 
   experim = experimC{i};
   disp(' ')
   disp(' --------------------------------------------------------')
   disp(['Experiment:        ' experim])

   switch Ropt
   case 1
     switch CS
     case 'R'
       CPS_R(experim)
     case 'Z'
       CPS_Z(experim,Gopt)
     case 'RZ'
       CPS_RZ(experim,Gopt)
     case {'XZ','YZ'}
       CPS_XZ(experim,Gopt)
     case 'XYZ'
       CPS_XYZ(experim,Gopt)
     end

   case 2
     switch CS
     case 'R'
       CVPM_R(experim)
     case 'Z'
       CVPM_Z(experim)
     case 'RZ'
       CVPM_RZ(experim)
     case {'XZ','YZ'}
       CVPM_XZ(experim)
     case 'XYZ'
       CVPM_XYZ(experim)
     end

   case 3
     switch CS
     case 'R'
       CPS_R( experim)
       CVPM_R(experim)
     case 'Z'
       CPS_Z( experim,Gopt)
       CVPM_Z(experim)
     case 'RZ'
       CPS_RZ( experim,Gopt)
       CVPM_RZ(experim)
     case {'XZ','YZ'}
       CPS_XZ( experim,Gopt)
       CVPM_XZ(experim)
     case 'XYZ'
       CPS_XYZ( experim,Gopt)
       CVPM_XYZ(experim)
     end
   end
 end
