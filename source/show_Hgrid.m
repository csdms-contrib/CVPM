 function [] = show_Hgrid(xfwL,xfeL,MtypL,xf,X,CS,C)

% Displays the horizontal (X,Y, or R) grid.
% ______________________________________________

% Notation:

%   CS = coordinate system
%   C  = the specific coordinate (X,Y or R).

 pos = set_screen2(0);

% Note: This function relies on Matlab plot functions.
% ______________________________________________

 N = length(X) - 1;

% define colors

 c = {[1 0 0],[0 0 1],[1 0.5 0],[0.5 1 0.8],[0 0.5 0.5],[0 0.5 0],[1 1 0],[0.95 0.95 0.95]};

% set limits

 Xmin = min(X);
 Xmax = max(X);
 jmin = find(xfwL <= Xmin,1,'last');
 jmax = find(xfeL >= Xmax,1,'first');

% show layers

 if ~strcmp(C,'Y')
   figure('position',pos)
   if ~strcmp(CS,'R')
     subplot(2,2,1)
   end
 else
   subplot(2,2,2)
 end

 for j=jmin:jmax
   switch fix(MtypL(j))
   case 1                       % test
     pc = c{1};
   case {2,3}                   % ice
     pc = c{2};
   case {4,5,6}                 % quartz, feldspar, mica-dominated igneouos/metamorphic rocks
     pc = c{3};
   case {7,8,9}                 % pyroxene, amphibole, olivine-dominated igneous/metamorphic rocks
     pc = c{4};
   case {10,11,12,13,14,15}     % sedimentary rocks
     pc = c{5};
   case {20,21,22,23,24,25}     % organic-rich materials
     pc = c{6};
   case {30,31,32,33,34,35}     % borehole fluids
     pc = c{7};
   case 99
     pc = c{8};
   otherwise
     disp(' ')
     disp('Error: show_Hgrid')
     disp('  patch color not defined for this MtypL')
     pause
   end
   pXX = [0 1 1 0];
   pX = [xfwL(j) xfwL(j) xfeL(j) xfeL(j)];
   patch(pX,pXX,pc,'Edgecolor','none')
   hold on
 end
 v = axis;
 xlabel([C ' (m) '])
 title('CV interfaces and grid points ')

 line([xfwL(1) xfwL(1)],v(3:4),'color','k','linewidth',4)
 for j=jmin:jmax
   line([xfeL(j) xfeL(j)],v(3:4),'color','k','linewidth',4)
 end

% show CV interfaces

 for i=1:N+1
   line([xf(i) xf(i)],v(3:4),'color',[0.8 0.8 0.8],'linewidth',1)
 end

% show CV grid pts

 gpts = 0.5*(v(4)-v(3))*ones(size(X));

 for i=1:N+1
   plot(X,gpts,'ko','MarkerFaceColor','k')
 end
 zoom on
