 function [] = show_HVgrid(xfwL,xfeL,xf,X,zfuL,zfdL,zf,Z,Mtyp,CS,C)

% Displays 2-D horizontal/vertical grid.

% Notation:

%   CS = coordinate system
%   C  = the specific coordinate (X,Y or R).

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

 N = length(X) - 1;
 M = length(Z) - 1;

% define colors

 c = {[1 0 0],[0 0 1],[1 0.5 0],[0.5 1 0.8],[0 0.5 0.5],[0 0.5 0],[1 1 0],[0.95 0.95 0.95]};

% set limits

 Xmin = min(X);
 Xmax = max(X);
 Zmin = min(Z);
 Zmax = max(Z);

 jmin = find(xfwL <= Xmin,1,'last');
 jmax = find(xfeL >= Xmax,1,'first');
 kmin = find(zfuL <= Zmin,1,'last');
 kmax = find(zfdL >= Zmax,1,'first');

% show layers

 figure('position',pos)
 set(gca,'YDir','reverse')
 hold on
 zoom on

 switch CS
 case 'RZ'
   v = [   0 Xmax Zmin Zmax];
 otherwise
   v = [Xmin Xmax Zmin Zmax];
 end
 axis(v)

 for j=1:N
   for k=1:M
     switch fix(Mtyp(k,j))
     case 1                         % test
       pc = c{1};
     case {2,3}                     % ice
       pc = c{2};
     case {4,5,6}                   % quartz, feldspar, mica-dominated igneous/metamorphic rocks
       pc = c{3};
     case {7,8,9}                   % pyroxene, amphibole, olive-dominated igneous/metamorphic rocks
       pc = c{4};
     case {10,11,12,13,14,15}       % sedimentary rocks
       pc = c{5};
     case {20,21,22,23,24,25}       % organic-rich materials
       pc = c{6};
     case {30,31,32,33,34,35}       % borehole fluids
       pc = c{7};
     otherwise
       disp(' ')
       disp('Error: show_Zgrid')
       disp('  patch color not defined for this Mtyp')
     end
     pX = [xf(j) xf(j+1) xf(j+1) xf(j)];
     pZ = [zf(k) zf(k) zf(k+1) zf(k+1)];
     patch(pX,pZ,pc,'Edgecolor','none')
   end
 end

 for j=jmin:jmax
   line([xfeL(j) xfeL(j)],v(3:4),'color','k','linewidth',4)
 end

 for k=kmin:kmax
   line([min(xfwL) max(xfeL)],[zfdL(k) zfdL(k)],'color','k','linewidth',4)
 end

% show CV interfaces

 for j=1:N+1
   line([xf(j) xf(j)],[Zmin Zmax],'color',[0.8 0.8 0.8],'linewidth',1)
 end

 for k=1:M+1
   line([Xmin Xmax],[zf(k) zf(k)],'color',[0.8 0.8 0.8],'linewidth',1)
 end

 xlabel([C ' (m) '])
 ylabel('Z (m) ')
 title('CV interfaces and grid points ')

% show CV grid pts

 for j=1:N+1
   switch j
   case 1
     x = xf(1);
   case N+1
     x = xf(N+1);
   otherwise
     x = xf(j) + 0.5*(xf(j+1) - xf(j));
   end
   gpts = x * ones(size(Z));
   plot(gpts,Z,'ko','MarkerFaceColor','k','markersize',4)
 end
