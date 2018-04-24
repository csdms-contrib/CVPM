% makeIC_Z.m

% Make an initial-condition file for the 1-D vertical CVPM case.

% This is currently setup to create a file containing the conditions just
% before a well is spudded.  It is assumed the vertical temperature profiles
% are the same everywhere prior to constructing the reserve pit and drilling
% pad.
% _______________________________________________________________

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
