% view_Z.m

% CVPM viewer for CVPM 1-D vertical simulations.

% Note: This function relies on Matlab plot functions.
% _____________________________________________________________

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
 
 dt_pause = 0.05;      % pause duration (sec)

 disp(' ')
 experim = input('Type name of experiment: ', 's');

% import fields created by CVPM 

 Ifile = ['CVPMout/' experim '_cvpm.mat'];

 load(Ifile)

 M     = length(Z) - 1;
 tgrid = tarr;
 Nout  = length(tgrid);
 Zmin  = min(Z);
 Zmax  = max(Z);
 secs_tunit = dtunit2secs(t_units);

 disp(' ')
 disp(['Zmax = ' num2str(Zmax)])
 disp(' ')
 zflag = input('Do you want to restrict the display to smaller z-values? [0 or 1]: ');
 if zflag
   Zmax = input('Type new Zmax: ');
   kp   = find(Z <= Zmax,1,'last');
   ZZ   = Z(1:kp);
   Zmax = max(ZZ);
 else
   kp   = M+1;
   ZZ   = Z;
 end

% > Display results

 figure('position',pos)
 for i=1:Nout

   switch i
   case Nout
     lw = 1.5;
   otherwise
     lw = 1;
   end

   T     = Tarr(    :,i);
   phi_i = phi_iarr(:,i);
   phi_u = phi_uarr(:,i);
   K     = Karr(    :,i);
   C     = Carr(    :,i);
   kappa = K ./ C;

   dTdz      = NaN*ones(size(T));
   dTdz(2:M) = (T(3:M+1)-T(1:M-1)) ./ (Z(3:M+1)-Z(1:M-1));

   if zflag
     T     = T(    1:kp);
     dTdz  = dTdz( 1:kp);
     phi_i = phi_i(1:kp);
     phi_u = phi_u(1:kp);
     K     = K(    1:kp);
     C     = C(    1:kp);
     kappa = kappa(1:kp);
   end

   ax(1) = subplot(2,3,1);
   plot(T,ZZ,'b','linewidth',lw)
   hold on 
   grid on
   zoom on
   set(gca,'YDir','reverse')
   xlabel('Temperature ($^\circ$C)','interpreter','latex')
   ylabel('Depth (m)','interpreter','latex')
   title(['Temperature~~(time = ' num2str(tgrid(i)) ' ' t_units ')'],'interpreter','latex')

   ax(2) = subplot(2,3,2);
   plot(1000*dTdz,ZZ,'r','linewidth',lw)
   hold on 
   grid on
   zoom on
   set(gca,'YDir','reverse')
   xlabel('$\partial T/\partial z$ (mK~m$^{-1}$)','interpreter','latex')
   title('Temperature Gradient','interpreter','latex')

   ax(3) = subplot(2,3,3);
   plot(K,ZZ,'color',[0 0.5 0],'linewidth',lw)
   hold on 
   grid on
   zoom on
   set(gca,'YDir','reverse')
   xlabel('$K$~(W~m$^{-1}$~K$^{-1}$)','interpreter','latex')
   title('Bulk Thermal Conductivity','interpreter','latex')

   ax(4) = subplot(2,3,4);
   plot(phi_i,ZZ,'b','linewidth',lw)
   hold on
   plot(phi_u,ZZ,'color',[0 0.5 0],'linewidth',lw)
   grid on
   zoom on
   set(gca,'YDir','reverse')
   xlabel('$\phi_i$ and $\phi_u$','interpreter','latex')
   ylabel('Depth (m)','interpreter','latex')
   title('Water Volume Fractions','interpreter','latex')
   v = axis;
   v(1) = 0;
   axis(v)
   if i == 1
     sleg = {' $\phi_i$ ',' $\phi_u$ '};
     legend(sleg,'location','southeast','interpreter','latex','fontsize',14)
   end

   ax(5) = subplot(2,3,5);
   plot(C/1e06,ZZ,'r','linewidth',lw)
   hold on 
   zoom on
   grid on
   set(gca,'YDir','reverse')
   xlabel('$C$ (MJ m$^{-3}$ K$^{-1}$) ','interpreter','latex')
   title('Volumetric Heat Capacity','interpreter','latex')

   ax(6) = subplot(2,3,6);
   plot(1e06*kappa,ZZ,'m','linewidth',lw)
   hold on 
   grid on
   zoom on
   set(gca,'YDir','reverse')
   xlabel('$10^6 \kappa$ (m$^2$ s$^{-1}$) ','interpreter','latex')
   title('Thermal Diffusivity','interpreter','latex')
   linkaxes(ax,'y')
   pause(dt_pause)
 end
 pause

% > Find errors for test cases

 if length(experim) < 4
   Exper = [experim '   '];
 else
   Exper = experim;
 end

 if strcmp(Exper(1:4),'Test')

   errorArr = NaN*ones(size(Tarr));

   Ts0 = Tsarr(1);

% find errors

   for i=1:Nout
     t     = secs_tunit * tgrid(i);
     T     = Tarr(:,i);
     Ts    = Tsarr(i);
     qb    = qbarr(i);
     C     = Carr(:,i);
     Tanal = Tz_analytic(experim,t,Ts0,Ts,qb,Z,Dz,dz,Km0,C,S0,hs);

     errorArr(:,i) = Tanal - T;
   end
   errorArr(M+1,:) = NaN;           % T(M+1) typically is just an estimate

% establish limits

   minQS    = min(min(QSarr));
   maxQS    = max(max(QSarr));
   minJ     = min(min(Jarr));
   maxJ     = max(max(Jarr));
   minJnet  = min(min(Jnetarr));
   maxJnet  = max(max(Jnetarr));
   minError = min(min(errorArr));
   maxError = max(max(errorArr));

   figure('position',pos)
   disp(' ')

   for i=1:Nout

     switch i
     case Nout
       lw = 1.5;
     otherwise
       lw = 1;
     end

     T      = Tarr(:,i);
     J      = Jarr(:,i);
     Jnet   = Jnetarr(:,i);
     error  = errorArr(:,i);
     maxErr = max(abs(error));
     J      = 1000*J;           % mW/m^2
     Jnet   = 1000*Jnet;        % mW/m^2
     QS     = 1000*QSarr(:,i);  % mW/m^2
 
     ax(1) = subplot(2,3,1);
     plot(T,Z,'r','linewidth',lw)
     hold on
     grid on
     zoom on
     set(gca,'YDir','reverse')
     v = axis;
     v(3:4) = [Zmin Zmax];
     axis(v)
     xlabel('$T$ ($^\circ$C)','interpreter','latex')
     ylabel('Depth (m)','interpreter','latex')
     title(['Temperature~~(time = ' num2str(tgrid(i)) ' ' t_units ')'],'interpreter','latex')

     ax(2) = subplot(2,3,2);
     plot(J,Z,'b','linewidth',lw)
     hold on
     grid on
     zoom on
     set(gca,'YDir','reverse')
     v(1)   = floor(1000*minJ);     % mW/m^2
     v(2)   = ceil( 1000*maxJ);     % mW/m^2
     v(3:4) = [Zmin Zmax];
     axis(v)
     xlabel('$J$ (mW~m$^{-2}$)','interpreter','latex')
     title('Heat Flux','interpreter','latex')

     ax(3) = subplot(2,3,3);
     plot(Jnet,Z,'color',[0 0.5 0],'linewidth',lw)
     hold on
     grid on
     zoom on
     set(gca,'YDir','reverse')
     v(1)   = floor(1000*minJnet);  % mW/m^2
     v(2)   = ceil( 1000*maxJnet);  % mW/m^2
     v(3:4) = [Zmin Zmax];
     if all(isfinite(v))
       axis(v)
     end
     xlabel('$J_\mathrm{net}$ (mW~m$^{-2}$)','interpreter','latex')
     title('Net Heat Flux','interpreter','latex')

     ax(4) = subplot(2,3,4);
     plot(error,Z,'m','linewidth',lw)
     hold on
     grid on
     zoom on
     set(gca,'YDir','reverse')
     x = floor(log10(maxError));
     v(1)   = 10^x * floor(10^(-x) * minError);
     v(2)   = 10^x * ceil( 10^(-x) * maxError);
     v(3:4) = [Zmin Zmax];
     if all(isfinite(v))
       axis(v)
     end
     line([0 0],v(3:4),'color','k','linestyle','--')
     xlabel('$\epsilon$ (K)','interpreter','latex')
     ylabel('Depth (m)','interpreter','latex')
     title('Error','interpreter','latex')

     ax(5) = subplot(2,3,5);
     plot(QS,Z,'r','linewidth',lw)
     hold on
     grid on
     zoom on
     set(gca,'YDir','reverse')
     v(1)   = floor(1000*minJnet);  % mW/m^2
     v(2)   = ceil( 1000*maxJnet);  % mW/m^2
     v(3:4) = [Zmin Zmax];
     if all(isfinite(v))
       axis(v)
     end
     xlabel('QS~(mW~m$^{-2}$)','interpreter','latex')
     title('Heat Production','interpreter','latex')

     ax(6) = subplot(2,3,6);
     plot(QS-Jnet,Z,'b','linewidth',lw)
     hold on
     grid on
     zoom on
     set(gca,'YDir','reverse')
     v(2)   = ceil( 1000*maxJnet)/5;
     v(1)   = -v(2);
     v(3:4) = [Zmin Zmax];
     if all(isfinite(v))
       axis(v)
     end
     xlabel('QS - $J_\mathrm{net}$~(mW~m$^{-2}$)','interpreter','latex')
     title('QS - $J_\mathrm{net}$','interpreter','latex')
     linkaxes(ax,'y')

     disp(['max(Error) = ' num2str(maxErr) ' K at t = ' num2str(tgrid(i)) ' '  t_units])
     pause(dt_pause)
   end
 end
