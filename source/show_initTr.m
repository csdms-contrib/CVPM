 function [] = show_initTr(rf,R,T,J,phi,K,C,rfwL,rfeL)

% Displays the initial temperature field (radial).

 pos = set_screen(0);

% Note: This function relies on Matlab plot functions.
% ______________________________________________

% Written by:

%   Gary Clow
%   Institute of Arctic and Alpine Research
%   University of Colorado
%   Boulder, Colorado USA
%   Email: gary.clow@colorado.edu
% ______________________________________________

 N = length(T) - 1;

% parameters

 lw = 1.0;      % line width for layers

% set limits

 Rmin = min(R);
 Rmax = max(R);
 jmin = find(rfwL <= Rmin,1,'last');
 jmax = find(rfeL >= Rmax,1,'first');

% show temperature profile

 figure('position',pos)
 subplot(2,3,1)
 plot(R,T,'r','linewidth',2)
 hold on
 grid on
 v = axis;
 v(1:2) = [Rmin Rmax];
 axis(v)
 xlabel('Radius (m)','interpreter','latex')
 ylabel('$T$ ($^\circ$C)','interpreter','latex')
 title('Initial Temperature Field','interpreter','latex')

% show temperature gradient

 dTdr = NaN*ones(size(T));

 for j=1:N
   dTdr(j) = (T(j+1)-T(j)) / (R(j+1)-R(j));
 end

 subplot(2,3,2)
 plot(R,dTdr,'b','linewidth',2)
 hold on
 grid on
 v = axis;
 v(1:2) = [Rmin Rmax];
 axis(v)
 xlabel('Radius (m)','interpreter','latex')
 ylabel('$\partial T / \partial r$ (mK m$^{-1}$)','interpreter','latex')
 title('Temperature Gradient','interpreter','latex')

% show diffusive heat-flux through the CV interfaces

 subplot(2,3,3)
 plot(rf,1000*J,'m','linewidth',2)
 hold on
 grid on
 v = axis;
 v(1:2) = [Rmin Rmax];
 axis(v)
 xlabel('Radius (m)','interpreter','latex')
 ylabel('$\tilde{J}$ (mW m$^{-1}$)','interpreter','latex')
 title('Integrated Heat Flux','interpreter','latex')

% show porosity

 subplot(2,3,4)
 plot(R,phi,'m','linewidth',2)
 hold on
 grid on
 v = axis;
 v(1:2) = [Rmin Rmax];
 axis(v)
 xlabel('Radius (m)','interpreter','latex')
 ylabel('$\phi$','interpreter','latex')
 title('Porosity','interpreter','latex')

% show bulk conductivity

 subplot(2,3,5)
 plot(R,K,'color',[0 0.5 0],'linewidth',2)
 hold on
 grid on
 v = axis;
 v(1:2) = [Rmin Rmax];
 axis(v)
 xlabel('Radius (m)','interpreter','latex')
 ylabel('$K$~(W~m$^{-1}$~K$^{-1}$)','interpreter','latex')
 title('Bulk Conductivity','interpreter','latex')

% show volumetric heat capacity

 subplot(2,3,6)
 plot(R,C/1e06,'color','r','linewidth',2)
 hold on
 grid on
 v = axis;
 v(1:2) = [Rmin Rmax];
 axis(v)
 xlabel('Radius (m)','interpreter','latex')
 ylabel('$C$ (MJ m$^{-3}$ K$^{-1}$) ','interpreter','latex')
 title('Volumetric Heat Capacity','interpreter','latex')

% show layers

 subplot(2,3,1)
 v = axis;
 line([rfwL(1) rfwL(1)],v(3:4),'color','k','linewidth',lw)
 for j=jmin:jmax
   line([rfeL(j) rfeL(j)],v(3:4),'color','k','linewidth',lw)
 end
 subplot(2,3,2)
 v = axis;
 line([rfwL(1) rfwL(1)],v(3:4),'color','k','linewidth',lw)
 for j=jmin:jmax
   line([rfeL(j) rfeL(j)],v(3:4),'color','k','linewidth',lw)
 end
 subplot(2,3,3)
 v = axis;
 line([rfwL(1) rfwL(1)],v(3:4),'color','k','linewidth',lw)
 for j=jmin:jmax
   line([rfeL(j) rfeL(j)],v(3:4),'color','k','linewidth',lw)
 end
 subplot(2,3,4)
 v = axis;
 line([rfwL(1) rfwL(1)],v(3:4),'color','k','linewidth',lw)
 for j=jmin:jmax
   line([rfeL(j) rfeL(j)],v(3:4),'color','k','linewidth',lw)
 end
 subplot(2,3,5)
 v = axis;
 line([rfwL(1) rfwL(1)],v(3:4),'color','k','linewidth',lw)
 for j=jmin:jmax
   line([rfeL(j) rfeL(j)],v(3:4),'color','k','linewidth',lw)
 end
 subplot(2,3,6)
 v = axis;
 line([rfwL(1) rfwL(1)],v(3:4),'color','k','linewidth',lw)
 for j=jmin:jmax
   line([rfeL(j) rfeL(j)],v(3:4),'color','k','linewidth',lw)
 end
