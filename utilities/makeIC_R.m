% makeIC_R.m

% Makes an initial-condition file for the 1-D radial CVPM case.

% This is typically done to simulate the drilling recovery phase of a 
% borehole.
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

% pickup the location of the working directory from the CVPM.config file

 [wdir,~,~,~] = input_config;
 cd(wdir)
% ______________________________________________

 disp(' ')
 disp('Type name of experiment from which to extract the temperature field:')
 experim = input('  ', 's');

% import fields created by CVPM

 Ifile = ['CVPMout/' experim '_cvpm.mat'];

 load(Ifile)

 secs_tunit = dtunit2secs(t_units);     % seconds per t_unit
 tt         = t / secs_tunit;

 disp(' ')
 disp(['Rmin = ' num2str(min(R)) ' (m)'])
 disp(['Rmax = ' num2str(max(R)) ' (m)'])
 disp(['t    = ' num2str(tt) ' (' t_units ')'])

 figure('position',pos)
 plot(R,T,'b','linewidth',2)
 hold on
 grid on
 zoom on
 xlabel('$R$ (m)','interpreter','latex')
 ylabel('$T (^\circ$C)','interpreter','latex')
 Sfile = convertS2latex(Ifile);
 title(['Final Temperature Field in ' Sfile ' at $t$ = ' num2str(tt) ' ' t_units],'interpreter','latex')

% assume temperatures inside borehole equal that at Rmin

 R = [0    R];
 T = [T(1) T];

 plot(R,T,'r','linewidth',1.5,'linestyle','--')

% create IC_file

 descript = ['final T-field from: ' experim];
 Imethod  = 'linear';
 IC       = T;
 varout   = 'descript R IC Imethod';

 fname    = [experim '_finalT'];
 Ofile    = ['ICs/' strtrim(fname) '.mat'];

 eval(['save ' Ofile ' ' varout]);
