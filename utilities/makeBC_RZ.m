% makeBC_RZ.m

% Makes a boundary-condition file for the 2-D cylindrical CVPM case.

% This is typically done to simulate the drilling phase of a borehole.
% ____________________________________________________________________

% > Set environment

 close all
 clear all
 format shortg
 colordef white
 pos = set_screen(3);

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

% parameters

 dz    = 2;         % resolution of vertical grid (m)
 dR    = 0.0092;    % drill pipe wall thickness
 Rdo   = 0.084;     % drill pipe outer radius
 Rh    = 0.155;     % open hole radius
 rho   = 1200;      % fluid density
 mf    = 36;        % fluid mass flux
 Te    = 30;        % fluid inlet temperatures
 T0    = -12;       % mean long-term surface temperature
 gamma = 0.04;      % mean geothermal gradient
 K     = 1.85;      % mean thermal conductivity of rocks
 kappa = 0.8e-06;   % mean thermal diffusivity of rocks

% inputs

 disp(' ')
 H = input('Type total borehole depth (m): ');

 disp(' ')
 tdays = input('Type number of days to drill hole: ');
 
% define depth grid

 Z = (dz:dz:H)';

% find annulus temperatures

 [td,Ta] = SB_drill_model(dR,Rdo,Rh,H,rho,mf,Te,T0,gamma,K,kappa,tdays,Z);

% find the drilling disturbance

 F = (T0 + gamma*Z);    % undisturbed formation profile

 dTa = zeros(size(Ta));

 for k=1:length(Z)
   dTa(k,:) = Ta(k,:) - F(k);
 end

 L      = isnan(dTa);
 dTa(L) = 0;            % fix dTa below the drill bit

% assume dTa at 0 m in the borehole is the same as first depth

 Z   = [0; Z];
 nZ  = length(Z);
 dTa = [dTa(1,:); dTa];

% assume dTa is zero when well is spudded (t = 0)

 td  = [0 td];
 dTa = [zeros(nZ,1) dTa];

% show drilling disturbance

 figure('position',pos)
 colormap jet
 contourf(td,Z,dTa)
 set(gca,'YDir','reverse')
 grid on
 xlabel('Time (days)')
 ylabel('Depth (m)')
 title('Drilling Disturbance','interpreter','latex')
 colorbar
 text(1.01,1.03,'$\Delta T_a$~(K)','units','normalized','interpreter','latex','fontsize',18)

% setup information needed by CVPM

 descript = ['drilling disturbance for a ' num2str(H) ' m borehole drilled in ' num2str(tdays) ' days'];
 t        = td';
 BC       = dTa';
 t_units  = 'days';
 BCtype   = 'T';
 Imethod  = 'linear';
 varout   = 'descript t t_units Z BCtype BC Imethod';

 disp(' ')
 fname = input('Type name of output file (without ext): ','s');

 Ofile = ['BCs/' fname '.mat'];

 eval(['save ' Ofile ' ' varout]);
