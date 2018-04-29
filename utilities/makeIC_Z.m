% makeIC_Z.m

% Make an initial-condition file for the 1-D vertical CVPM case.

% This is currently setup to create a file containing the conditions just
% before a well is spudded.  It is assumed the vertical temperature profiles
% are the same everywhere prior to constructing the reserve pit and drilling
% pad.
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

% > Set environment

 close all
 clear all
 format shortg
 colordef white
 pos = set_screen(0);

% set working directory

 [wdir,~,~,~] = input_config;
 cd(wdir)
% _______________________________

 disp(' ')
 disp('Type name of 1-D CVPM experiment used to simulate tundra temps')
 disp(' just prior to well construction:')
 experim = input('  ', 's');

% import fields created by CVPM

 Ifile = ['CVPMout/' experim '_cvpm.mat'];

 load(Ifile)

 figure('position',pos)
 plot(T,Z)
 grid on
 zoom on
 set(gca,'Ydir','reverse')
 xlabel('Temperature ($^\circ$C)','interpreter','latex')
 ylabel('Depth (m)','interpreter','latex')
 title('Temperature Field Just Prior to Spudding','interpreter','latex')
 pause

% create IC_file

 descript = [experim ', just prior to spudding the hole'];
 Imethod  = 'linear';
 IC       = T;
 varout   = 'descript Z IC Imethod';

 Ofile = ['ICs/IC_' experim '_z.mat'];

 eval(['save ' Ofile ' ' varout]);
