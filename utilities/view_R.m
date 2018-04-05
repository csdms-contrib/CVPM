% view_R.m

% CVPM viewer for CVPM 1-D radial simulations.

% Note: This function relies on Matlab plot functions.
% _____________________________________________________________

% > Set environment

 close all
 clear all
 format shortg
 colordef white
 pos = set_screen2(0);

% pickup the location of the working directory from the CVPM.config file

 [wdir,~,~,~] = input_config;
 addpath([wdir '/source'],'-end')
 addpath([wdir '/utilities'],'-end')
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

 N     = length(R) - 1;
 tgrid = tarr;
 Nout  = length(tgrid);
 Rmin  = min(R);
 Rmax  = max(R);
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
   kp   = N+1;
   RR   = R;
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

   dTdr = NaN*ones(size(T));
   for k=2:N
     dTdr(k) = (T(k+1)-T(k-1)) / (R(k+1)-R(k-1));
   end

   if rflag
     T     = T(    1:kp);
     dTdr  = dTdr( 1:kp);
     phi_i = phi_i(1:kp);
     phi_u = phi_u(1:kp);
     K     = K(    1:kp);
     C     = C(    1:kp);
     kappa = kappa(1:kp);
   end

   ax(1) = subplot(2,3,1);
   plot(RR,T,'b','linewidth',lw)
   hold on 
   grid on
   zoom on
   xlabel('Radius (m)','interpreter','latex')
   ylabel('$T$ ($^\circ$C)','interpreter','latex')
   title(['Temperature~~(time = ' num2str(tgrid(i)) ' ' t_units ')'],'interpreter','latex')

   ax(2) = subplot(2,3,2);
   plot(RR,dTdr,'r','linewidth',lw)
   hold on 
   grid on
   xlabel('Radius (m)','interpreter','latex')
   ylabel('$\partial T/\partial z$ (K~m$^{-1}$)','interpreter','latex')
   title('Temperature Gradient','interpreter','latex')

   ax(3) = subplot(2,3,3);
   plot(RR,K,'color',[0 0.5 0],'linewidth',lw)
   hold on 
   grid on
   xlabel('Radius (m)','interpreter','latex')
   ylabel('$K$~(W~m$^{-1}$~K$^{-1}$)','interpreter','latex')
   title('Bulk Thermal Conductivity','interpreter','latex')

   ax(4) = subplot(2,3,4);
   plot(RR,phi_i,'b','linewidth',lw)
   hold on
   grid on
   xlabel('Radius (m)','interpreter','latex')
   ylabel('$\phi_i$','Interpreter','latex')
   title('Volume Fraction of Ice','interpreter','latex')

   ax(5) = subplot(2,3,5);
   plot(RR,C/1e06,'r','linewidth',lw)
   hold on 
   grid on
   xlabel('Radius (m)','interpreter','latex')
   ylabel('$C$ (MJ m$^{-3}$ K$^{-1}$) ','interpreter','latex')
   title('Volumetric Heat Capacity','interpreter','latex')

   ax(6) = subplot(2,3,6);
   plot(RR,1e06*kappa,'m','linewidth',lw)
   hold on 
   grid on
   xlabel('Radius (m)','interpreter','latex')
   ylabel('$10^6 \kappa$ (m$^2$ s$^{-1}$) ','interpreter','latex')
   title('Thermal Diffusivity','interpreter','latex')
   linkaxes(ax,'x')
   pause(1)
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

   Ta0 = Taarr(1);

   for i=1:Nout

%   find errors

     t     = secs_tunit * tgrid(i);
     T     = Tarr(:,i);
     Ta    = Taarr(i);
     qo    = qoarr(i);
     C     = Carr(:,i);
     Tanal = Tr_analytic(experim,t,Ta0,Ta,qo,R,Dr,dr,Km0,C,S0,hs);

     errorArr(:,i) = Tanal - T;
   end
   errorArr(N+1,:) = NaN;           % T(N+1) typically is just an estimate

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
     plot(T,R,'r','linewidth',lw)
     hold on
     grid on
     set(gca,'YDir','reverse')
     v = axis;
     v(3:4) = [Rmin Rmax];
     axis(v)
     xlabel('$T$ ($^\circ$C)','interpreter','latex')
     ylabel('Depth (m)','interpreter','latex')
     title(['Temperature~~(time = ' num2str(tgrid(i)) ' ' t_units ')'],'interpreter','latex')

     ax(2) = subplot(2,3,2);
     plot(J,R,'b','linewidth',lw)
     hold on
     grid on
     set(gca,'YDir','reverse')
     v(1)   = floor(1000*minJ);
     v(2)   = ceil( 1000*maxJ);
     v(3:4) = [Rmin Rmax];
     axis(v)
     xlabel('$J$ (mW~m$^{-2}$)','interpreter','latex')
     title('Heat Flux','interpreter','latex')

     ax(3) = subplot(2,3,3);
     plot(Jnet,R,'color',[0 0.5 0],'linewidth',lw)
     hold on
     grid on
     set(gca,'YDir','reverse')
     v(1)   = floor(1000*minJnet);
     v(2)   = ceil( 1000*maxJnet);
     v(3:4) = [Rmin Rmax];
     if all(isfinite(v))
       axis(v)
     end
     xlabel('$J_\mathrm{net}$ (mW~m$^{-2}$)','interpreter','latex')
     title('Net Heat Flux','interpreter','latex')

     ax(4) = subplot(2,3,4);
     plot(error,R,'m','linewidth',lw)
     hold on
     grid on
     set(gca,'YDir','reverse')
     x = floor(log10(maxError));
     v(1)   = 10^x * floor(10^(-x) * minError);
     v(2)   = 10^x * ceil( 10^(-x) * maxError);
     v(3:4) = [Rmin Rmax];
     if all(isfinite(v))
       axis(v)
     end
     line([0 0],v(3:4),'color','k','linestyle','--')
     xlabel('$\epsilon$ (K)','interpreter','latex')
     ylabel('Depth (m)','interpreter','latex')
     title('Error','interpreter','latex')

     ax(5) = subplot(2,3,5);
     plot(QS,R,'r','linewidth',lw)
     hold on
     grid on
     set(gca,'YDir','reverse')
     v(1)   = floor(1000*minJnet);
     v(2)   = ceil( 1000*maxJnet);
     v(3:4) = [Rmin Rmax];
     if all(isfinite(v))
       axis(v)
     end
     xlabel('QS~(mW~m$^{-2}$)','interpreter','latex')
     title('Heat Production','interpreter','latex')

     ax(6) = subplot(2,3,6);
     plot(QS-Jnet,R,'b','linewidth',lw)
     hold on
     grid on
     set(gca,'YDir','reverse')
     v(2)   = ceil( 1000*maxJnet)/5;
     v(1)   = -v(2);
     v(3:4) = [Rmin Rmax];
     if all(isfinite(v))
       axis(v)
     end
     xlabel('QS - $J_\mathrm{net}$~(mW~m$^{-2}$)','interpreter','latex')
     title('QS - $J_\mathrm{net}$','interpreter','latex')
     linkaxes(ax,'y')

     disp(['max(Error) = ' num2str(maxErr) ' K at t = ' num2str(tgrid(i)) ' '  t_units])
     pause(1)
   end
 end
