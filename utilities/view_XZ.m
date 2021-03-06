% view_XZ.m

% CVPM viewer for CVPM 2-D cartesian simulations.
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

% Note: This function relies on Matlab plot functions.
% _____________________________________________________________

% Set environment

 close all
 clear all
 format shortg
 colordef white
 pos = set_screen2(0);

% plot options

 bflag = 0;     % bflag=1 shows vertical line at x=0, y=0
 pflag = 1;     % pflag=1 shows dT/dz, pflag=2 shows (T - Tinit)

% pickup the location of the working directory from the CVPM.config file

 [wdir,~,~,~] = input_config;
 addpath([wdir '/source'],'-end')
 addpath([wdir '/utilities'],'-end')
 cd(wdir)
% ______________________________________________
 
 disp(' ')
 experim = input('Type name of experiment: ', 's');

% import fields created by CVPM

 Ifile = ['CVPMout/' experim '_cvpm.mat'];

 load(Ifile)

 M     = length(Z) - 1;
 N     = length(X) - 1;
 tgrid = tarr;
 Nout  = length(tgrid);
 Zmin  = min(Z);
 Zmax  = max(Z);
 Xmin  = min(X);
 Xmax  = max(X);
 Tinit = Tarr(:,:,1);
 secs_tunit = dtunit2secs(t_units);

 disp(' ')
 disp(['Zmax = ' num2str(Zmax)])
 disp(' ')
 zflag = input('Do you want to restrict the display to shallower depths? [0 or 1]: ');
 if zflag
   Zmax  = input('Type new Zmax: ');
   kp    = find(Z <= Zmax,1,'last');
   ZZ    = Z(1:kp);
   Zmax  = max(ZZ);
   Tinit = Tinit(1:kp,:);
 else
   kp    = M+1;
   ZZ    = Z;
 end

% > Display results

 figure('position',pos)

 v = [Xmin Xmax Zmin Zmax];

 for i=1:Nout
   T     = Tarr(    :,:,i);
   phi_i = phi_iarr(:,:,i);
   phi_u = phi_uarr(:,:,i);
   K     = Karr(    :,:,i);
   C     = Carr(    :,:,i);

   dTdz = NaN*ones(size(T));
   dTdz(2:M,:) = (T(3:M+1,:)-T(1:M-1,:)) ./ repmat(Z(3:M+1)-Z(1:M-1),[1,N+1]);

   if zflag
     T     = T(    1:kp,:);
     dTdz  = dTdz( 1:kp,:);
     phi_i = phi_i(1:kp,:);
     phi_u = phi_u(1:kp,:);
     K     = K(    1:kp,:);
     C     = C(    1:kp,:);
   end

   kappa = K ./ C;
   delT  = T - Tinit;

   subplot(2,3,1)
   colormap jet
   contourf(X,ZZ,T)
   grid on
   set(gca,'Ydir','reverse')
   axis(v)
   xlabel('Horizontal Distance (m)','interpreter','latex')
   ylabel('Depth (m)','interpreter','latex')
   title(['Temperature Field  (' num2str(tgrid(i)) ' ' t_units ')'],'interpreter','latex')
   if bflag
     line([0 0],v(3:4),'color','w','linewidth',2)
   end
   colorbar

   subplot(2,3,2)
   colormap jet
   switch pflag
   case 1
     contourf(X,ZZ,1000*dTdz)
     title('Gradient, $\partial T/\partial z$ (mK~m$^{-1}$)','interpreter','latex')
   case 2
     contourf(X,ZZ,delT)
     title('$\Delta T$ (K)','interpreter','latex')
   end
   grid on
   set(gca,'YDir','reverse')
   axis(v)
   xlabel('Horizontal Distance (m)','interpreter','latex')
   if bflag
     line([0 0],v(3:4),'color','w','linewidth',2)
   end
   colorbar

   subplot(2,3,3)
   colormap jet
   contourf(X,ZZ,K)
   grid on
   set(gca,'YDir','reverse')
   axis(v)
   xlabel('Horizontal Distance (m)','interpreter','latex')
   title('Conductivity, $K$~(W~m$^{-1}$~K$^{-1}$)','interpreter','latex')
   colorbar

   subplot(2,3,4)
   colormap jet
   contourf(X,ZZ,phi_i)
   grid on
   set(gca,'YDir','reverse')
   axis(v)
   xlabel('Horizontal Distance (m)','interpreter','latex')
   ylabel('Depth (m)','interpreter','latex')
   title('Volume Fraction of Ice, $\phi_i$','interpreter','latex')
   colorbar

   subplot(2,3,5)
   colormap jet
   contourf(X,ZZ,log10(C/1e06))
   grid on
   set(gca,'Ydir','reverse')
   axis(v)
   xlabel('Horizontal Distance (m)','interpreter','latex')
   title('Heat Capacity, log[$C$] (MJ~m$^{-3}$~K$^{-1}$)','interpreter','latex')
   colorbar

   subplot(2,3,6)
   colormap jet
   contourf(X,ZZ,1e06*kappa)
   grid on
   set(gca,'Ydir','reverse')
   axis(v)
   xlabel('Horizontal Distance (m)','interpreter','latex')
   title('Diffusivity, $10^6 \kappa$ (m$^2$ s$^{-1}$) ','interpreter','latex')
   colorbar
   pause(0.01)
 end
 pause

% > 3-D plots

 dt5  = round((max(tgrid) - min(tgrid)) / 5);
 tvec = min(tgrid) + (0:1:5)*dt5;

 figure('position',pos)
 colormap jet
 T = Tarr;
 if zflag
   T = T(1:kp,:,:);
 end
 Tp = permute(T,[3 2 1]);
 [X3,t3,Z3] = meshgrid(X,tgrid,ZZ);
 xslice = [];
 tslice = tvec;
 zslice = [];
 h = slice(X3,t3,Z3,Tp,xslice,tslice,zslice,'spline');
 set(h,'edgecolor','none')
 set(gca,'Zdir','reverse')
 xlabel('Horizontal Distance (m)','interpreter','latex')
 ylabel(t_units,'interpreter','latex')
 zlabel('Depth (m)','interpreter','latex')
 title('Temperature ($^\circ$C)','interpreter','latex')
 colorbar
 view(165,50)
 rotate3d on
 pause

 figure('position',pos)
 cmap = colormap(flipud(jet));
 nc   = size(cmap,1);
 cmap(1,:) = [1 50/256 1];
 colormap(cmap)
 phi_i = phi_iarr;
 if zflag
   phi_i = phi_i(1:kp,:,:);
 end
 phi_ip = permute(phi_i,[3 2 1]);
 [X3,t3,Z3] = meshgrid(X,tgrid,ZZ);
 xslice = [];
 tslice = tvec;
 zslice = [];
 h = slice(X3,t3,Z3,phi_ip,xslice,tslice,zslice);
 set(h,'edgecolor','none')
 set(gca,'Zdir','reverse')
 xlabel('Horizontal Distance (m)','interpreter','latex')
 ylabel(t_units,'interpreter','latex')
 zlabel('Depth (m)','interpreter','latex')
 title('Volume Fraction of Ice, $\phi_i$','interpreter','latex')
 colorbar
 view(165,50)
 rotate3d on
 pause

 figure('position',pos)
 colormap jet
 phi_u = phi_uarr;
 if zflag
   phi_u = phi_u(1:kp,:,:);
 end
 phi_up = permute(phi_u,[3 2 1]);
 [X3,t3,Z3] = meshgrid(X,tgrid,ZZ);
 xslice = [];
 tslice = tvec;
 zslice = [];
 h = slice(X3,t3,Z3,phi_up,xslice,tslice,zslice);
 set(h,'edgecolor','none')
 set(gca,'Zdir','reverse')
 xlabel('Horizontal Distance (m)','interpreter','latex')
 ylabel(t_units,'interpreter','latex')
 zlabel('Depth (m)','interpreter','latex')
 title('Volume Fraction of Unfrozen Water, $\phi_\ell$','interpreter','latex')
 colorbar
 view(165,50)
 rotate3d on
 pause

 figure('position',pos)
 colormap jet
 K = Karr;
 if zflag
   K = K(1:kp,:,:);
 end
 Kp = permute(K,[3 2 1]);
 [X3,t3,Z3] = meshgrid(X,tgrid,ZZ);
 xslice = [];
 tslice = tvec;
 zslice = [];
 h = slice(X3,t3,Z3,Kp,xslice,tslice,zslice);
 set(h,'edgecolor','none')
 set(gca,'Zdir','reverse')
 xlabel('Horizontal Distance (m)','interpreter','latex')
 ylabel(t_units,'interpreter','latex')
 zlabel('Depth (m)','interpreter','latex')
 title('Conductivity, $k$~(W~m$^{-1}$~K$^{-1}$)','interpreter','latex')
 colorbar
 view(165,50)
 rotate3d on
 pause

 figure('position',pos)
 colormap jet
 C = Carr;
 if zflag
   C = C(1:kp,:,:);
 end
 Cp = permute(C,[3 2 1]);
 [X3,t3,Z3] = meshgrid(X,tgrid,ZZ);
 xslice = [];
 tslice = tvec;
 zslice = [];
 h = slice(X3,t3,Z3,log10(Cp/1e06),xslice,tslice,zslice);
 set(h,'edgecolor','none')
 set(gca,'Zdir','reverse')
 xlabel('Horizontal Distance (m)','interpreter','latex')
 ylabel(t_units,'interpreter','latex')
 zlabel('Depth (m)','interpreter','latex')
 title('Diffusivity, $10^6 \kappa$ (m$^2$ s$^{-1}$) ','interpreter','latex')
 colorbar
 view(165,50)
 rotate3d on
 pause

 figure('position',pos)
 colormap jet
 K = Karr;
 C = Carr;
 if zflag
   K = K(1:kp,:,:);
   C = C(1:kp,:,:);
 end
 kappa  = K ./ C;
 kappap = permute(kappa,[3 2 1]);
 [X3,t3,Z3] = meshgrid(X,tgrid,ZZ);
 xslice = [];
 tslice = tvec;
 zslice = [];
 h = slice(X3,t3,Z3,kappap*1e06,xslice,tslice,zslice);
 set(h,'edgecolor','none')
 set(gca,'Zdir','reverse')
 xlabel('Horizontal Distance (m)','interpreter','latex')
 ylabel(t_units,'interpreter','latex')
 zlabel('Depth (m)','interpreter','latex')
 title('Diffusivity, $10^6 \kappa$ (m$^2$ s$^{-1}$) ','interpreter','latex')
 colorbar
 view(165,50)
 rotate3d on
 pause


% > Find errors for test cases

 if length(experim) < 4
   Exper = [experim '   '];
 else
   Exper = experim;
 end

 if strcmp(Exper(1:4),'Test')

% find errors

   errorArr = NaN*ones(size(Tarr));
   Ts0      = Tsarr(1);

   for i=1:Nout
     t  = secs_tunit * tgrid(i);
     Ts = Tsarr( :,i);
     qb = qbarr( :,i);
     T  = Tarr(:,:,i);
     C  = Carr(:,:,i);

     Tanal = NaN*ones(M+1,N+1);
     for j=1:N+1
       Tanal(:,j) = Tz_analytic(experim,t,Ts0,Ts(j),qb(j),Z,Dz,dz,Km0(:,j),C(:,j),S0(:,j),hs(:,j));
     end
     errorArr(:,:,i) = Tanal - T;
   end

   errorArr(M+1,:,:) = NaN;         % set boundaries to NaN
   errorArr(:,1,  :) = NaN;
   errorArr(:,N+1,:) = NaN;

% show errors

   figure('position',pos)
   disp(' ')

   for i=1:Nout
     T     = Tarr(    :,:,i);
     error = errorArr(:,:,i);
     maxErr = max(max(abs(error)));

     subplot(2,2,1)
     colormap jet
     contourf(X,Z,T)
     grid on
     set(gca,'Ydir','reverse')
     axis(v)
     xlabel('Horizontal Distance (m)','interpreter','latex')
     ylabel('Depth (m)','interpreter','latex')
     title(['Temperature Field  (' num2str(tgrid(i)) ' ' t_units ')'],'interpreter','latex')
     colorbar

     subplot(2,2,3)
     colormap jet
     contourf(X,Z,abs(error))
     grid on
     set(gca,'Ydir','reverse')
     axis(v)
     xlabel('Horizontal Distance (m)','interpreter','latex')
     ylabel('Depth (m)','interpreter','latex')
     title('$|$Error$|$ (K)','interpreter','latex')
     colorbar

     disp(['max(Error) = ' num2str(maxErr) ' K at t = ' num2str(tgrid(i)) ' '  t_units])
     pause(0.2)
   end
 end
