% makeIC_RZ.m

% Makes an initial-condition file for the 2-D cylindrical CVPM case.

% This is typically done to simulate the drilling recovery phase of a 
% borehole.
% _______________________________________________________________

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

% Written by:

%   Gary Clow
%   Institute of Arctic and Alpine Research
%   University of Colorado
%   Boulder, Colorado USA
%   Email: gary.clow@colorado.edu
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
 colormap jet
 contourf(R,Z,T)
 grid on
 zoom on
 set(gca,'Ydir','reverse')
 xlabel('Radial Distance (m)','interpreter','latex')
 ylabel('Depth (m)','interpreter','latex')
 Sfile = convertS2latex(Ifile);
 title(['Final Temperature Field in ' Sfile ' at $t$ = ' num2str(tt) ' ' t_units],'interpreter','latex')
 colorbar
 pause

% assume temperatures inside borehole equals that at Rmin

 R = [0      R];
 T = [T(:,1) T];

 figure('position',pos)
 colormap jet
 contourf(R,Z,T)
 grid on
 zoom on
 set(gca,'Ydir','reverse')
 v = axis;
 line([R(2) R(2)],v(3:4),'color','w')
 xlabel('Radial Distance (m)','interpreter','latex')
 ylabel('Depth (m)','interpreter','latex')
 title(['Final Temperature Field in ' Sfile ' at $t$ = ' num2str(tt) ' ' t_units],'interpreter','latex')
 colorbar

% create IC_file

 descript = ['final T-field from: ' experim];
 Imethod  = 'linear';
 IC       = T;
 varout   = 'descript R Z IC Imethod';

 fname    = [experim '_finalT'];
 Ofile    = ['ICs/' strtrim(fname) '.mat'];

 eval(['save ' Ofile ' ' varout]);
