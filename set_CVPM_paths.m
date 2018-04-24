% set_CVPM_paths.m

% Sets the paths to required CVPM directories.
% ____________________________________________

% > Set CVPM working directory

% *** change the following statement to satisfy your working environment ***

 wdir = '~/thermal/numer/CVPM_v1.1';

% ************************************

% > Add paths

 addpath( wdir,                '-end')
 addpath([wdir '/AK_projects'],'-end')
 addpath([wdir '/source'],     '-end')
 addpath([wdir '/utilities'],  '-end')
 savepath
