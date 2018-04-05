 function [] = show_initTz(zf,Z,T,J,phi,K,C,zfuL,zfdL)

% Displays the initial temperature field (vertical).

 pos = set_screen2(0);

% Note: This function relies on Matlab plot functions.
% ______________________________________________

% Written by:

%   Gary Clow
%   Institute of Arctic and Alpine Research
%   University of Colorado
%   Boulder, Colorado USA
%   Email: gary.clow@colorado.edu
% ______________________________________________

 M = length(T) - 1;

% parameters

 lw = 1.0;      % line width for layers

% set limits

 Zmin = min(Z);
 Zmax = max(Z);
 jmin = find(zfuL <= Zmin,1,'last');
 jmax = find(zfdL >= Zmax,1,'first');

% show temperature profile

 figure('position',pos)
 subplot(2,3,1)
 plot(T,Z,'r','linewidth',2)
 hold on
 grid on
 set(gca,'YDir','reverse')
 v = axis;
 v(3:4) = [Zmin Zmax];
 axis(v)
 xlabel('Temperature ($^\circ$C)','interpreter','latex')
 ylabel('Depth (m)','interpreter','latex')
 title('Initial Temperature Field','interpreter','latex')

% show temperature gradient

 dTdz = NaN*ones(size(T));

 for j=1:M
   dTdz(j) = (T(j+1)-T(j)) / (Z(j+1)-Z(j));
 end

 subplot(2,3,2)
 plot(1000*dTdz,Z,'b','linewidth',2)
 hold on
 grid on
 v = axis;
 v(3:4) = [Zmin Zmax];
 axis(v)
 set(gca,'YDir','reverse')
 xlabel('$\partial T / \partial z$ (mK m$^{-1}$)','interpreter','latex')
 title('Temperature Gradient','interpreter','latex')

% show diffusive heat-flux through the CV interfaces

 subplot(2,3,3)
 plot(1000*J,zf,'m','linewidth',2)
 hold on
 grid on
 v = axis;
 v(3:4) = [Zmin Zmax];
 axis(v)
 set(gca,'YDir','reverse')
 xlabel('$J$ (mW m$^{-2}$)','interpreter','latex')
 title('Diffusive Heat Flux','interpreter','latex')

% show porosity

 subplot(2,3,4)
 plot(phi,Z,'m','linewidth',2)
 hold on
 grid on
 set(gca,'YDir','reverse')
 v = axis;
 v(3:4) = [Zmin Zmax];
 axis(v)
 xlabel('$\phi$','interpreter','latex')
 ylabel('Depth (m)','interpreter','latex')
 title('Porosity','interpreter','latex')

% show bulk conductivity

 subplot(2,3,5)
 plot(K,Z,'color',[0 0.5 0],'linewidth',2)
 hold on
 grid on
 set(gca,'YDir','reverse')
 v = axis;
 v(3:4) = [Zmin Zmax];
 axis(v)
 xlabel('$K$~(W~m$^{-1}$~K$^{-1}$)','interpreter','latex')
 title('Bulk Conductivity','interpreter','latex')

% show volumetric heat capacity

 subplot(2,3,6)
 plot(C/1e06,Z,'color','r','linewidth',2)
 hold on
 grid on
 set(gca,'YDir','reverse')
 v = axis;
 v(3:4) = [Zmin Zmax];
 axis(v)
 xlabel('$C$ (MJ m$^{-3}$ K$^{-1}$) ','interpreter','latex')
 title('Volumetric Heat Capacity','interpreter','latex')

% show layers

 subplot(2,3,1)
 v = axis;
 line(v(1:2),[zfuL(1) zfuL(1)],'color','k','linewidth',lw)
 for j=jmin:jmax
   line(v(1:2),[zfdL(j) zfdL(j)],'color','k','linewidth',lw)
 end
 subplot(2,3,2)
 v = axis;
 line(v(1:2),[zfuL(1) zfuL(1)],'color','k','linewidth',lw)
 for j=jmin:jmax
   line(v(1:2),[zfdL(j) zfdL(j)],'color','k','linewidth',lw)
 end
 subplot(2,3,3)
 v = axis;
 line(v(1:2),[zfuL(1) zfuL(1)],'color','k','linewidth',lw)
 for j=jmin:jmax
   line(v(1:2),[zfdL(j) zfdL(j)],'color','k','linewidth',lw)
 end
 subplot(2,3,4)
 v = axis;
 line(v(1:2),[zfuL(1) zfuL(1)],'color','k','linewidth',lw)
 for j=jmin:jmax
   line(v(1:2),[zfdL(j) zfdL(j)],'color','k','linewidth',lw)
 end
 subplot(2,3,5)
 v = axis;
 line(v(1:2),[zfuL(1) zfuL(1)],'color','k','linewidth',lw)
 for j=jmin:jmax
   line(v(1:2),[zfdL(j) zfdL(j)],'color','k','linewidth',lw)
 end
 subplot(2,3,6)
 v = axis;
 line(v(1:2),[zfuL(1) zfuL(1)],'color','k','linewidth',lw)
 for j=jmin:jmax
   line(v(1:2),[zfdL(j) zfdL(j)],'color','k','linewidth',lw)
 end
