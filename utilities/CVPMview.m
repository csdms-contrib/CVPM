% CVPMview.m

% View the results from a CVPM simulation.
% ______________________________________________

% Written by:

%   Gary Clow
%   Institute of Arctic and Alpine Research
%   University of Colorado
%   Boulder, Colorado USA
%   Email: gary.clow@colorado.edu
% ______________________________________________

% pickup the location of the working directory from the CVPM.config file

 [wdir,~,~,~] = input_config;
 cd(wdir)
% __________________________________________________________________________

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
