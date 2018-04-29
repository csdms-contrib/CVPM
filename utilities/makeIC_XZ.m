% makeIC_XZ.m

% Makes an initial-condition file for the 2-D cartesian CVPM case by 
% extracting the final temperature field from a previous simulation.  If
% the previous simulation was a 1-D vertical experiment, the extracted
% temperature field is replicated across the 2-D problem domain.  If the 
% previous simulation was a 2-D cartesian experiment, the new field is the
% same as the final field of that experiment.
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
 Sfile = convertS2latex(Ifile);

 secs_tunit = dtunit2secs(t_units);     % seconds per t_unit
 tt         = t / secs_tunit;

 disp(' ')
 disp(['CS   = ' CS])
 if strcmp(CS,'XZ')
   disp(['Xmin = ' num2str(min(X)) ' (m)'])
   disp(['Xmax = ' num2str(max(X)) ' (m)'])
 end
 disp(['Zmin = ' num2str(min(Z)) ' (m)'])
 disp(['Zmax = ' num2str(max(Z)) ' (m)'])
 disp(['t    = ' num2str(tt) ' (' t_units ')'])

 switch CS
 case 'Z'
   X  = [-5000 0 5000];
   N  = length(X);
   Tz = T;
   T  = repmat(Tz,1,N);	% assume initial temp field is independent of X

 case 'XZ'
   figure('position',pos)
   colormap jet
   contourf(X,Z,T)
   grid on
   zoom on
   set(gca,'Ydir','reverse')
   xlabel('Horizontal Distance (m)','interpreter','latex')
   ylabel('Depth (m)','interpreter','latex')
   title(['Final Temperature Field in ' Sfile ' at $t$ = ' num2str(tt) ' ' t_units],'interpreter','latex')
   colorbar
   pause
 end

 figure('position',pos)
 colormap jet
 contourf(X,Z,T)
 grid on
 zoom on
 set(gca,'Ydir','reverse')
 xlabel('Horizontal Distance (m)','interpreter','latex')
 ylabel('Depth (m)','interpreter','latex')
 title(['Final Temperature Field in ' Sfile ' at $t$ = ' num2str(tt) ' ' t_units],'interpreter','latex')
 colorbar

% create IC_file

 descript = ['final T-field from: ' experim];
 Imethod  = 'linear';
 IC       = T;
 varout   = 'descript X Z IC Imethod';

 fname    = [experim '_finalT'];
 Ofile    = ['ICs/' strtrim(fname) '.mat'];

 eval(['save ' Ofile ' ' varout]);
