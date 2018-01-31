% makeBC_RZ.m

% Makes a boundary-condition file for the 2-D cylindrical CVPM case.

% This is typically done to simulate the drilling/conditioning phases of a
% borehole.

% Notes: 

%   (1) At present, this code is only setup to create a drilling disturbance
%       BC based on output from AK_make_varphi_p1td.
%   (2) Times are "local" time, ie. t=0 when drill reaches this depth. 
% ____________________________________________________________________

% > Set environment

 close all
 clear all
 format shortg
 colordef white

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
 well_code = input('Type well_code: ','s');

% load drilling disturbance varphi (part 1) file

 fn = ['AK_projects/drilling_disturb/varphi_' well_code '_p1td_cell'];
 load(fn)

 nZ = size(varphiCell,1);
 td = varphiCell{1,8};      % time since well was spudded (global time) <all depths use the same time array>
 nt = length(td);

% preallocate arrays

 Z      = NaN*ones(nZ,1);
 DeltaT = NaN*ones(nZ,nt);

 for k=1:nZ               % depth loop ------------->

   Z(k)      = varphiCell{k,1};
   delTstar  = varphiCell{k,6};     % scaling temperature
   varphi_td = varphiCell{k,9};     % varphi on global timescale

% find warming at borehole wall

   DeltaT(k,:) = delTstar * varphi_td;
 end

% assume DeltaT at 0 m in the borehole is the same as that at 1 m depth

 Z      = [0; Z];
 nZ     = length(Z);
 DeltaT = [DeltaT(1,:); DeltaT];

% assume DeltaT is zero when well is spudded (t = 0)

 td     = [0 td];
 DeltaT = [zeros(nZ,1) DeltaT];

% setup information needed by CVPM

 descript = [well_code  ', drilling disturbance, part 1'];
 t        = td';
 BC       = DeltaT';
 t_units  = 'days';
 BCtype   = 'T';
 Imethod  = 'linear';
 varout   = 'descript t t_units Z BCtype BC Imethod';

 disp(' ')
 fname = input('Type name of output file (without ext): ','s');

 Ofile = ['BCs/' fname '.mat'];

 eval(['save ' Ofile ' ' varout]);
