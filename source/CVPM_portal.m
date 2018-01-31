% CVPM_portal.m

% This is the portal to the CVPM permafrost modeling system.

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

% Written by:

%   Gary Clow
%   Institute of Arctic and Alpine Research
%   University of Colorado
%   Boulder, Colorado USA
%   Email: gary.clow@colorado.edu
% ______________________________________________

 [wdir,CS,opt,experimC] = input_config;

 disp(' ')
 disp(['Working Directory: ' wdir])
 disp(' ')
 disp(['Coordinate System: ' CS])

% set working directory

 cd(wdir)

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
