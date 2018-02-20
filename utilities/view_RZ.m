% view_RZ.m

% CVPM viewer for CVPM 2-D cylindrical simulations.

% Note: This function relies on Matlab plot functions.
% _____________________________________________________________

% > Set environment

 close all
 clear all
 format shortg
 colordef white
 pos = set_screen(0);

% plot mode

 pflag = 2;     % pflag=1 show dT/dr, pflag=2 shows (T - Tinit)

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
 experim = input('Type name of experiment: ', 's');

% import fields created by CVPM

 Ifile = ['CVPMout/' experim '_cvpm.mat'];

 load(Ifile)

 upperBC_type = BCtypes{1};
 lowerBC_type = BCtypes{2};
 innerBC_type = BCtypes{3};
 outerBC_type = BCtypes{4};

 N     = length(R) - 1;
 M     = length(Z) - 1;
 tgrid = tarr;
 Nout  = length(tgrid);
 Rmin  = min(R);
 switch outerBC_type
 case 'T'
   Rmax = R(N+1);
 case 'q'
   Rmax = R(N);
 end
 Zmin  = min(Z);
 switch lowerBC_type
 case 'T'
   Zmax = Z(M+1);
 case 'q'
   Zmax = Z(M);
 end
 Tinit = Tarr(:,:,1);
 secs_tunit = dtunit2secs(t_units);

 disp(' ')
 disp(['Rmax = ' num2str(Rmax)])
 disp(' ')
 rflag = input('Do you want to restrict the display to smaller r-values? [0 or 1]: ');
 if rflag
   Rmax = input('Type new Rmax: ');
   kp   = find(R <= Rmax,1,'last');
   RR   = R(1:kp);
   Rmax = max(RR);
 else
   kp = N+1;
   RR = R;
 end
 Tfar = Tinit(:,N);      % initial far-field values
 Tfar = repmat(Tfar,1,kp);

% > Display results

 figure('position',pos)

 v = [0 Rmax Zmin Zmax];

 for i=1:Nout

   T     = Tarr(    :,:,i);
   phi_i = phi_iarr(:,:,i);
   phi_u = phi_uarr(:,:,i);
   K     = Karr(    :,:,i);
   C     = Carr(    :,:,i);

   dTdr = NaN*ones(size(T));
   for k=2:N
     dTdr(:,2:N) = (T(:,3:N+1)-T(:,1:N-1)) ./ repmat(R(3:N+1)-R(1:N-1),[M+1 1]);
   end

   if rflag
     T     = T(    :,1:kp);
     dTdr  = dTdr( :,1:kp);
     phi_i = phi_i(:,1:kp);
     phi_u = phi_u(:,1:kp);
     K     = K(    :,1:kp);
     C     = C(    :,1:kp);
   end

   kappa = K ./ C;
   delT  = T - Tfar;

   ax(1) = subplot(2,3,1);
   colormap jet
   contourf(RR,Z,T)
   grid on
   set(gca,'Ydir','reverse')
   axis(v)
   xlabel('Radial Distance (m)','interpreter','latex')
   ylabel('Depth (m)','interpreter','latex')
   title(['Temperature Field  (' num2str(tgrid(i)) ' ' t_units ')'],'interpreter','latex')
   colorbar

   ax(2) = subplot(2,3,2);
   colormap jet
   switch pflag
   case 1
     contourf(RR,Z,dTdr)
     title('$\partial T/\partial r$ (K~m$^{-1}$)','interpreter','latex')
   case 2
     contourf(RR,Z,delT)
     title('$\Delta T$ (K)','interpreter','latex')
   end
   grid on
   set(gca,'YDir','reverse')
   axis(v)
   xlabel('Radial Distance (m)','interpreter','latex')
   colorbar

   ax(3) = subplot(2,3,3);
   colormap jet
   contourf(RR,Z,K)
   grid on
   set(gca,'YDir','reverse')
   axis(v)
   xlabel('Radial Distance (m)','interpreter','latex')
   title('Conductivity, $K$~(W~m$^{-1}$~K$^{-1}$)','interpreter','latex')
   colorbar

   ax(4) = subplot(2,3,4);
   colormap jet
   contourf(RR,Z,phi_i)
   grid on
   set(gca,'YDir','reverse')
   axis(v)
   xlabel('Radial Distance (m)','interpreter','latex')
   ylabel('Depth (m)','interpreter','latex')
   title('Volume Fraction of Ice, $\phi_i$','interpreter','latex')
   colorbar

   ax(5) = subplot(2,3,5);
   colormap jet
   contourf(RR,Z,C/1e06)
   grid on
   set(gca,'Ydir','reverse')
   axis(v)
   xlabel('Radial Distance (m)','interpreter','latex')
   title('Heat Capacity, $C$ (MJ/m$^3$ K)','interpreter','latex')
   colorbar

   ax(6) = subplot(2,3,6);
   colormap jet
   contourf(RR,Z,1e06*kappa)
   grid on
   set(gca,'Ydir','reverse')
   axis(v)
   xlabel('Radial Distance (m)','interpreter','latex')
   title('Diffusivity, $10^6 \kappa$ (m$^2$ s$^{-1}$) ','interpreter','latex')
   colorbar
   linkaxes(ax,'xy')
   pause(0.01)
 end
 pause

 if length(experim) < 4
   Exper = [experim '   '];
 else
   Exper = experim;
 end

 if strcmp(Exper(1:5),'Test_')

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
     contourf(R,Z,T)
     grid on
     set(gca,'Ydir','reverse')
     axis(v)
     xlabel('Radial Distance (m)','interpreter','latex')
     ylabel('Depth (m)','interpreter','latex')
     title(['Temperature Field  (' num2str(tgrid(i)) ' ' t_units ')'],'interpreter','latex')
     colorbar

     subplot(2,2,3)
     colormap jet
     contourf(R,Z,abs(error))
     grid on
     set(gca,'Ydir','reverse')
     axis(v)
     xlabel('Radial Distance (m)','interpreter','latex')
     ylabel('Depth (m)','interpreter','latex')
     title('$|$Error$|$ (K)','interpreter','latex')
     colorbar

     disp(['max(Error) = ' num2str(maxErr) ' K at t = ' num2str(tgrid(i)) ' '  t_units])
     pause(0.2)
   end
   pause
 end

% > Show temperature profile at min(R)

 figure('position',pos)

 for i=1:Nout
   f    = (i-1)/(Nout-1);
   s{i} = [f 0 1-f];
 end

 Tfar = Tinit(:,N);      % initial far-field values

 for i=1:Nout
   T    = Tarr(:,1,i);
   delT = T - Tfar;

   ax(1) = subplot(1,2,1);
   plot(T,Z,'color',s{i},'linewidth',1.5)
   hold on
   grid on
   set(gca,'YDir','reverse')
   xlabel('Temperature ($^\circ$C)','interpreter','latex')
   ylabel('Depth (m)','interpreter','latex')
   title(['Time = ' num2str(tgrid(i)) ' ' t_units ',  r = ' num2str(100*min(R)) ' cm'])

   ax(2) = subplot(1,2,2);
   plot(delT,Z,'color',s{i},'linewidth',1.5)
   hold on
   grid on
   set(gca,'YDir','reverse')
   xlabel('$\Delta T$ (K)','interpreter','latex')
   linkaxes(ax,'y')
   pause(0.01)
 end
 pause

% > Show final temperature profiles at several r-values

 figure('position',pos)

 dj   = max(1,floor(N/20));
 jmax = round(0.7*N);
 ic   = 1;
 Tfar = Tarr(:,N,Nout);

 for j=1:dj:jmax
   f        = (j-1)/(jmax-1);
   lc       = [f 0 1-f];
   sleg{ic} = [' $r = $' num2str(R(j)) ' m'];

   T    = Tarr(:,j,Nout);
   delT = T - Tfar;

   ax(1) = subplot(1,2,1);
   plot(T,Z,'color',lc,'linewidth',1.5)
   hold on
   grid on
   set(gca,'YDir','reverse')
   xlabel('Temperature ($^\circ$C)','interpreter','latex')
   ylabel('Depth (m)','interpreter','latex')
   title(['Temperature Profiles at ' num2str(tgrid(Nout)) ' ' t_units])

   ax(2) = subplot(1,2,2);
   semilogx(delT,Z,'color',lc,'linewidth',1.5)
   hold on
   grid on
   set(gca,'YDir','reverse')
   xlabel('$\Delta T$ (K)','interpreter','latex')
   title('Multiple r-values')
   linkaxes(ax,'y')
   pause(0.01)
   ic = ic + 1;
 end
 subplot(1,2,1)
 v    = axis;
 v(4) = max(Z);
 axis(v)
 subplot(1,2,2)
 v    = axis;
 v(4) = max(Z);
 axis(v)
 legend(sleg,'location','southwest','interpreter','latex','fontsize',14)
