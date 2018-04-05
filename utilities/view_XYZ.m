% view_XYZ.m

% CVPM viewer for CVPM 3-D cartesian simulations.

% Note: This function relies on Matlab plot functions.
% _____________________________________________________________

% > Set environment

 close all
 clear all
 format shortg
 colordef white
 pos  = set_screen2(0);
 pos2 = pos;
 pos2(3) = 1.2*pos2(3);

% plot options

 bflag = 1;     % bflag=1 shows vertical line at x=0, y=0
 pflag = 1;     % pflag=1 show dT/dz, pflag=2 shows (T - Tinit)

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

 M     = length(Z) - 1;
 N     = length(X) - 1;
 L     = length(Y) - 1;
 tgrid = tarr;
 Nout  = length(tgrid);
 Zmin  = min(Z);
 Zmax  = max(Z);
 Xmin  = min(X);
 Xmax  = max(X);
 Ymin  = min(Y);
 Ymax  = max(Y);
 Tinit = Tarr(:,:,:,1);
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
   Tinit = Tinit(1:kp,:,:);
 else
   kp    = M+1;
   ZZ    = Z;
 end

% select where to show cross-sections

 if Xmin ~= 0 && Ymin ~= 0       % this is probably a drill pad experiment
   junk = abs(X);
   xp   = find(junk == min(junk),1,'first');
   junk = abs(Y);
   yp   = find(junk == min(junk),1,'first');
 else
   midX = (Xmin + Xmax)/2;
   junk = abs(X - midX);
   xp   = find(junk == min(junk),1,'first');
   midY = (Ymin + Ymax)/2;
   junk = abs(Y - midY);
   yp   = find(junk == min(junk),1,'first');
 end

% > Display an XZ cross-section

 disp(' ')
 disp(['Displaying results for a vertical slice at y = ' num2str(Y(yp))])

 figure('position',pos)

 vx = [Xmin Xmax Zmin Zmax];
 vy = [Ymin Ymax Zmin Zmax];

 for i=1:Nout
   T     = Tarr(    :,:,yp,i);
   phi_i = phi_iarr(:,:,yp,i);
   phi_u = phi_uarr(:,:,yp,i);
   K     = Karr(    :,:,yp,i);
   C     = Carr(    :,:,yp,i);

   dTdz = NaN*ones(size(T));
   dTdz(2:M,:) = (T(3:M+1,:)-T(1:M-1,:)) ./ repmat(Z(3:M+1)-Z(1:M-1),[1 N+1]);

   if zflag
     T     = T(    1:kp,:);
     dTdz  = dTdz( 1:kp,:);
     phi_i = phi_i(1:kp,:);
     phi_u = phi_u(1:kp,:);
     K     = K(    1:kp,:);
     C     = C(    1:kp,:);
   end

   kappa = K ./ C;
   delT  = T - Tinit(:,:,yp);

   subplot(2,3,1)
   colormap jet
   contourf(X,ZZ,T)
   grid on
   set(gca,'Ydir','reverse')
   axis(vx)
   xlabel('$X$ (m)','interpreter','latex')
   ylabel('Depth (m)','interpreter','latex')
   title(['Temperature Field  (' num2str(tgrid(i)) ' ' t_units ')'],'interpreter','latex')
   if bflag
     line([0 0],vx(3:4),'color','w','linewidth',2)
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
   axis(vx)
   xlabel('$X$ (m)','interpreter','latex')
   if bflag
     line([0 0],vx(3:4),'color','w','linewidth',2)
   end
   colorbar

   subplot(2,3,3)
   colormap jet
   contourf(X,ZZ,K)
   grid on
   set(gca,'YDir','reverse')
   axis(vx)
   xlabel('$X$ (m)','interpreter','latex')
   title('Conductivity, $K$~(W~m$^{-1}$~K$^{-1}$)','interpreter','latex')
   colorbar

   subplot(2,3,4)
   colormap jet
   contourf(X,ZZ,phi_i)
   grid on
   set(gca,'YDir','reverse')
   axis(vx)
   xlabel('$X$ (m)','interpreter','latex')
   ylabel('Depth (m)','interpreter','latex')
   title('Volume Fraction of Ice, $\phi_i$','interpreter','latex')
   colorbar

   subplot(2,3,5)
   colormap jet
   contourf(X,ZZ,C/1e06)
   grid on
   set(gca,'Ydir','reverse')
   v = axis;
   v(2:4) = [Xmax Zmin Zmax];
   axis(vx)
   xlabel('$X$ (m)','interpreter','latex')
   title('Heat Capacity, $C$ (MJ/m$^3$ K)','interpreter','latex')
   colorbar

   subplot(2,3,6)
   colormap jet
   contourf(X,ZZ,1e06*kappa)
   grid on
   set(gca,'Ydir','reverse')
   axis(vx)
   xlabel('$X$ (m)','interpreter','latex')
   title('Diffusivity, $10^6 \kappa$ (m$^2$ s$^{-1}$) ','interpreter','latex')
   colorbar
   pause(0.01)
 end
 pause

% > Display a YZ cross-section

 disp(' ')
 disp(['Displaying results for a vertical slice at x = ' num2str(X(xp))])

 figure('position',pos)

 v = [Ymin Ymax Zmin Zmax];

 for i=1:Nout
   T     = squeeze(Tarr(    :,xp,:,i));
   phi_i = squeeze(phi_iarr(:,xp,:,i));
   phi_u = squeeze(phi_uarr(:,xp,:,i));
   K     = squeeze(Karr(    :,xp,:,i));
   C     = squeeze(Carr(    :,xp,:,i));

   dTdz = NaN*ones(size(T));
   dTdz(2:M,:) = (T(3:M+1,:)-T(1:M-1,:)) ./ repmat(Z(3:M+1)-Z(1:M-1),[1 L+1]);

   if zflag
     T     = T(    1:kp,:);
     dTdz  = dTdz( 1:kp,:);
     phi_i = phi_i(1:kp,:);
     phi_u = phi_u(1:kp,:);
     K     = K(    1:kp,:);
     C     = C(    1:kp,:);
   end

   kappa = K ./ C;
   delT  = T - squeeze(Tinit(:,xp,:));

   subplot(2,3,1)
   colormap jet
   contourf(Y,ZZ,T)
   grid on
   set(gca,'Ydir','reverse')
   axis(vy)
   xlabel('$Y$ (m)','interpreter','latex')
   ylabel('Depth (m)','interpreter','latex')
   title(['Temperature Field  (' num2str(tgrid(i)) ' ' t_units ')'],'interpreter','latex')
   if bflag
     line([0 0],vy(3:4),'color','w','linewidth',2)
   end
   colorbar

   subplot(2,3,2)
   colormap jet
   switch pflag
   case 1
     contourf(Y,ZZ,1000*dTdz)
     title('Gradient, $\partial T/\partial z$ (mK~m$^{-1}$)','interpreter','latex')
   case 2
     contourf(Y,ZZ,delT)
     title('$\Delta T$ (K)','interpreter','latex')
   end
   grid on
   set(gca,'YDir','reverse')
   axis(vy)
   xlabel('$Y$ (m)','interpreter','latex')
   if bflag
     line([0 0],vy(3:4),'color','w','linewidth',2)
   end
   colorbar

   subplot(2,3,3)
   colormap jet
   contourf(Y,ZZ,K)
   grid on
   set(gca,'YDir','reverse')
   axis(vy)
   xlabel('$Y$ (m)','interpreter','latex')
   title('Conductivity, $K$~(W~m$^{-1}$~K$^{-1}$)','interpreter','latex')
   colorbar

   subplot(2,3,4)
   colormap jet
   contourf(Y,ZZ,phi_i)
   grid on
   set(gca,'YDir','reverse')
   axis(vy)
   xlabel('$Y$ (m)','interpreter','latex')
   ylabel('Depth (m)','interpreter','latex')
   title('Volume Fraction of Ice, $\phi_i$','interpreter','latex')
   colorbar

   subplot(2,3,5)
   colormap jet
   contourf(Y,ZZ,C/1e06)
   grid on
   set(gca,'Ydir','reverse')
   axis(vy)
   xlabel('$Y$ (m)','interpreter','latex')
   title('Heat Capacity, $C$ (MJ/m$^3$ K)','interpreter','latex')
   colorbar

   subplot(2,3,6)
   colormap jet
   contourf(Y,ZZ,1e06*kappa)
   grid on
   set(gca,'Ydir','reverse')
   axis(vy)
   xlabel('$Y$ (m)','interpreter','latex')
   title('Diffusivity, $10^6 \kappa$ (m$^2$ s$^{-1}$) ','interpreter','latex')
   colorbar
   pause(0.01)
 end
 pause

% > Dislay slice planes

% Notes:
%   (1) Slice requires the dimensions of the array T be permuted in such a
%       way that it basically goes (Ty,Tx,Tz).
%   (2) The axes displayed by the following code are essentially index 
%       numbers, not actual values.  Thus a value of 10 on the x-axis
%       corresponds to the 10th element of X, NOT X = 10 m.  Haven't 
%       figured out how to fix this yet.

 dx1 = abs(X -   0);
 dx2 = abs(X - 110);
 dxm = abs(X - Xmax);
 dy1 = abs(Y -   0);
 dy2 = abs(Y -  50);
 dym = abs(Y - Ymax);
 dz1 = abs(Z -   0);
 dz2 = abs(Z -  20);
 dz3 = abs(Z -  50);
 dzm = abs(Z - Zmax);

 ix1  = find(dx1 == min(dx1),1,'first');
 ix2  = find(dx2 == min(dx2),1,'first');
 ixm  = find(dxm == min(dxm),1,'first');
 iy1  = find(dy1 == min(dy1),1,'first');
 iy2  = find(dy2 == min(dy2),1,'first');
 iym  = find(dym == min(dym),1,'first');
 iz1  = find(dz1 == min(dz1),1,'first');
 iz2  = find(dz2 == min(dz2),1,'first');
 iz3  = find(dz3 == min(dz3),1,'first');
 izm  = find(dzm == min(dzm),1,'first');

 sx = [ixm];
 sy = [iym];
 sz = [iz1,iz2,iz3,izm];

 sx2 = [ixm];
 sy2 = [iy1,iym];
 sz2 = [izm];

 sx3 = [ix1,ixm];
 sy3 = [iym];
 sz3 = [izm];

 figure('position',pos2)
 colormap jet
 for i=1:Nout

   T = squeeze(Tarr(:,:,:,i));
   if zflag
     T = T(1:kp,:,:);
   end

   Tp = permute(T,[3 2 1]);

   subplot(1,3,1)
   h1 = slice(Tp,sx,sy,sz);
   set(h1,'FaceColor','interp','EdgeColor','none')
   set(gca,'ZDir','reverse')
   xlabel('X ')
   ylabel('Y ')
   zlabel('Z ')
   title(['Temperature Field  (' num2str(tgrid(i)) ' ' t_units ')  '])
   axis tight
   box on
   view(-38.5,16)
   camzoom(1.2)
   colorbar('SouthOutside')
%   camproj perspective

   subplot(1,3,2)
   h2  = slice(Tp,sx2,sy2,sz2);
   set(h2,'FaceColor','interp','EdgeColor','none')
   set(gca,'ZDir','reverse')
   xlabel('X ')
   ylabel('Y ')
   axis tight
   box on
   view(-38.5,16)
   camzoom(1.2)
   colorbar('SouthOutside')
%   camproj perspective

   subplot(1,3,3)
   h2  = slice(Tp,sx3,sy3,sz3);
   set(h2,'FaceColor','interp','EdgeColor','none')
   set(gca,'ZDir','reverse')
   xlabel('X ')
   ylabel('Y ')
   axis tight
   box on
   view(-38.5,16)
   camzoom(1.2)
   colorbar('SouthOutside')
%   camproj perspective
   pause(0.2)
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
   Ts0      = Tsarr(1,1);

   for i=1:Nout
     t  = secs_tunit * tgrid(i);
     Ts = Tsarr( :,:,i);
     qb = qbarr( :,:,i);
     T  = Tarr(:,:,:,i);
     C  = Carr(:,:,:,i);

     Tanal = NaN*ones(M+1,N+1,L+1);
     for jy=1:L+1
       for j=1:N+1
         Tanal(:,j,jy) = Tz_analytic(experim,t,Ts0,Ts(j,jy),qb(j,jy),Z,Dz,dz,Km0(:,j,jy),C(:,j,jy),S0(:,j,jy),hs(:,j,jy));
       end
     end
     errorArr(:,:,:,i) = Tanal - T;
   end

   errorArr(M+1,:,:,  :) = NaN;         % set boundaries to NaN
   errorArr(:,1,  :,  :) = NaN;
   errorArr(:,N+1,:,  :) = NaN;
   errorArr(:,:  ,1,  :) = NaN;
   errorArr(:,:  ,L+1,:) = NaN;   

   figure('position',pos)
   disp(' ')

   for i=1:Nout
     T        = Tarr(    :,:,:,i);
     error    = abs(errorArr(:,:,:,i));
     maxErr   = max(max(max(error)));
     maxXYerr = squeeze(max(error,[],1));   % max error along each column

     subplot(2,2,1)
     colormap jet
     contourf(X,Z,T(:,:,yp))
     grid on
     set(gca,'Ydir','reverse')
     axis(vx)
     xlabel('$X$ (m)','interpreter','latex')
     ylabel('Depth (m)','interpreter','latex')
     title(['Temperature Field  (' num2str(tgrid(i)) ' ' t_units ')'],'interpreter','latex')
     colorbar

     subplot(2,2,2)
     colormap jet
     contourf(X,Y,maxXYerr)
     xlabel('$X$ (m)','interpreter','latex')
     ylabel('$Y$ (m)','interpreter','latex')
     title('max$|$Error$|$ (K)','interpreter','latex')
     colorbar

     subplot(2,2,3)
     colormap jet
     contourf(X,Z,error(:,:,yp))
     grid on
     set(gca,'Ydir','reverse')
     axis(vx)
     xlabel('$X$ (m)','interpreter','latex')
     ylabel('Depth (m)','interpreter','latex')
     title('$|$Error$|$ (K)','interpreter','latex')
     colorbar

     subplot(2,2,4)
     colormap jet
     contourf(Y,Z,squeeze(error(:,xp,:)))
     grid on
     set(gca,'Ydir','reverse')
     axis(vy)
     xlabel('$Y$ (m)','interpreter','latex')
     ylabel('Depth (m)','interpreter','latex')
     title('$|$Error$|$ (K)','interpreter','latex')
     colorbar

     disp(['max(Error) = ' num2str(maxErr) ' K at t = ' num2str(tgrid(i)) ' '  t_units])
     pause(0.2)
   end
   pause
 end

% > Show temperature profile at (xp,yp)

 disp(['Temperature profile at (x,y) = ' num2str(X(xp)) ',' num2str(Y(yp)) ' m'])

 if Nout >= 10
   dN = 2;
 else
   dN = 1;
 end

 for i=1:Nout
   f = (i-1)/(Nout-1);
   s{i} = [f 0 1-f];
 end

 nn     = length(t_units);
 t_unit = t_units(1:nn-1);
 ic     = 1;

 figure('position',pos)

 for i=1:dN:Nout
   T    = Tarr(1:kp,xp,yp,i);
   delT = T - Tinit(:,xp,yp);

   legstr{ic} = ['\hspace{0.2em} year ' num2str(tgrid(i)) '\hspace{0.3em}'];

   subplot(1,2,1)
   plot(T,ZZ,'color',s{i},'linewidth',1.5)
   hold on
   grid on
   set(gca,'YDir','reverse')
   xlabel('Temperature ($^\circ$C)','interpreter','latex')
   ylabel('Depth (m)','interpreter','latex')
   title(['Time = ' num2str(tgrid(i)) ' ' t_units])

   subplot(1,2,2)
   plot(delT,ZZ,'color',s{i},'linewidth',1.5)
   hold on
   grid on
   set(gca,'YDir','reverse')
   xlabel('$\Delta T$ (K)','interpreter','latex')
   pause(1)
   ic = ic + 1;
 end

 subplot(1,2,1)
 [hleg,hobj] = legend(legstr,'location','northeast','fontsize',14,'interpreter','latex');
 set(hleg,'linewidth',0.5,'edgecolor',[0.5 0.5 0.5])
