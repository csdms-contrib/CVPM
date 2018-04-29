 function [] = show_Vgrid(zfuL,zfdL,MtypL,zf,Z,CS)

% Displays the vertical grid.
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

 pos  = set_screen2(3);
 pos1 = [pos(1)-320 pos(2)+400 0.8*pos(3:4)];

% Note: This function relies on Matlab plot functions.
% ______________________________________________

 M = length(Z) - 1;

% define colors

 c = {[1 0 0],[0 0 1],[1 0.5 0],[0.5 1 0.8],[0 0.5 0.5],[0 0.5 0],[1 1 0],[0.95 0.95 0.95],[1 0 0]};

% set limits

 Zmin = min(Z);
 Zmax = max(Z);
 kmin = find(zfuL <= Zmin,1,'last');
 kmax = find(zfdL >= Zmax,1,'first');

% show layers

 switch CS
 case 'Z'
   figure('position',pos1)
 otherwise
   subplot(2,2,3)
 end

 for k=kmin:kmax
   switch fix(MtypL(k))
   case 1                       % test
     pc = c{1};
   case {2,3}                   % ice
     pc = c{2};
   case {4,5,6}                 % quartz, feldspar, mica-dominated igneous/metamorphic rocks
     pc = c{3};
   case {7,8,9}                 % pyroxene, amphibole, olivine-dominated igneous/metamorphic rocks
     pc = c{4};
   case {10,11,12,13,14,15}     % sedimentary rocks
     pc = c{5};
   case {20,21,22,23,24,25}     % organic-rich materials
     pc = c{6};
   case {30,31,32,33,34,35}     % borehole fluids
     pc = c{7};
   case {40,41,42,43,44}        % metals
     pc = c{8};
   otherwise
     disp(' ')
     disp('Error: show_Zgrid')
     disp('  patch color not defined for this MtypL')
   end
   pX = [0 1 1 0];
   pZ = [zfuL(k) zfuL(k) zfdL(k) zfdL(k)];
   patch(pX,pZ,pc,'Edgecolor','none')
   hold on
 end
 set(gca,'YDir','reverse')
 v = axis;
 ylabel('Z (m) ')
 title('CV interfaces and grid points ')

 line(v(1:2),[zfuL(1) zfuL(1)],'color','k','linewidth',4)
 for k=kmin:kmax
   line(v(1:2),[zfdL(k) zfdL(k)],'color','k','linewidth',4)
 end

% show CV interfaces

 for i=1:M+1
   line(v(1:2),[zf(i) zf(i)],'color',[0.8 0.8 0.8],'linewidth',1)
 end

% show CV grid pts

 gpts = 0.5*(v(2)-v(1))*ones(size(Z));

 for i=1:M+1
   plot(gpts,Z,'ko','MarkerFaceColor','k')
 end
 zoom on
