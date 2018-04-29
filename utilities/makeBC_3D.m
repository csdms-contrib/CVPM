% makeBC_3D.m

% Makes a boundary-condition file for the 3-D cartesian CVPM case.
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

% > Set environment

 close all
 clear all
 format shortg
 colordef white

% pickup the location of the working directory from the CVPM.config file

 [wdir,~,~,~] = input_config;
 cd(wdir)
% ______________________________________________

 disp(' ')
 disp('Available boundaries: ')
 disp('[1] XY face at min(Z), aka upper Z')
 disp('[2] XY face at max(Z), aka lower Z')
 disp('[3] YZ face at min(X)')
 disp('[4] YZ face at max(X)')
 disp('[5] XZ face at min(Y)')
 disp('[6] XZ face at max(Y)')
 Copt = input('Select boundary (1,2,3...): ');

 disp(' ')
 switch Copt
 case {1,2}
   X = input('Type X vector: ');
   Y = input('Type Y vector: ');
   N = length(X);
   M = length(Y);
 case {3,4}
   Y = input('Type Y vector: ');
   Z = input('Type Z vector: ');
   N = length(Y);
   M = length(Z);
 case {5,6}
   X = input('Type X vector: ');
   Z = input('Type Z vector: ');
   N = length(X);
   M = length(Z);
 end

 disp(' ')
 t  = input('Type time vector: ');
 nt = length(t);

 BC = NaN*ones(nt,M,N);

 disp(' ')
 t_units = input('Type time units: ','s');

 disp(' ')
 disp('Available boundary conditions: ')
 disp('[T] prescribed temperature')
 disp('[q] prescribed heat flux')
 BCtype = input('Select BC (T or q): ','s');

 disp(' ')
 disp('Available methods for creating the BC at each time-slice: ')
 disp('[1] uniform BC across the face.')
 disp('[2] one differing patch within a uniform background field.')
 disp('[3] two differing patches within a uniform background field.')
 Mopt = input('Select method (1,2,3...): ');

 switch Copt
 case {1,2}
   switch Mopt
   case 2
     disp(' ')
     cLL = input('Type X and Y locations of lower-left corner (m): ');
     cUR = input('Type X and Y locations of upper-right corner (m): ');
     x0  = cLL(1);
     y0  = cLL(2);
     xm  = cUR(1);
     ym  = cUR(2);
     Lx  = X >= x0 & X <= xm;
     Ly  = Y >= y0 & Y <= ym;
     Ly  = Ly';
     Lxy = logical(repmat(Lx,[M 1]) .* repmat(Ly,[1 N]));
   case 3
     disp(' ')
     cLL = input('Patch 1: type X and Y locations of lower-left corner (m): ');
     cUR = input('Patch 1: type X and Y locations of upper-right corner (m): ');
     x0  = cLL(1);
     y0  = cLL(2);
     xm  = cUR(1);
     ym  = cUR(2);
     Lx  = X >= x0 & X <= xm;
     Ly  = Y >= y0 & Y <= ym;
     Ly  = Ly';
     Lxy1 = logical(repmat(Lx,[M 1]) .* repmat(Ly,[1 N]));
     cLL = input('Patch 2: type X and Y locations of lower-left corner (m): ');
     cUR = input('Patch 2: type X and Y locations of upper-right corner (m): ');
     x0  = cLL(1);
     y0  = cLL(2);
     xm  = cUR(1);
     ym  = cUR(2);
     Lx  = X >= x0 & X <= xm;
     Ly  = Y >= y0 & Y <= ym;
     Ly  = Ly';
     Lxy2 = logical(repmat(Lx,[M 1]) .* repmat(Ly,[1 N]));
   end
 end

 disp(' ')
 for i=1:nt
   switch Mopt
   case 1
     val = input(['Type value at time ' num2str(t(i)) ' ' t_units ': ']);
     BC(i,:,:) = val;
   case 2
     valB = input(['Type background value at time ' num2str(t(i)) ' ' t_units ': ']);
     valP = input(['Type patch      value at time ' num2str(t(i)) ' ' t_units ': ']);
     A      = valB * ones(size(Lxy));
     A(Lxy) = valP
     BC(i,:,:) = permute(A,[3 1 2]);
   case 3
     valB  = input(['Type background value at time ' num2str(t(i)) ' ' t_units ': ']);
     valP1 = input(['Type patch 1    value at time ' num2str(t(i)) ' ' t_units ': ']);
     valP2 = input(['Type patch 2    value at time ' num2str(t(i)) ' ' t_units ': ']);
     A       = valB * ones(size(Lxy1));
     A(Lxy1) = valP1;
     A(Lxy2) = valP2
     BC(i,:,:) = permute(A,[3 1 2]);
   end
 end

 disp(' ')
 Imethod = input('Type CVPM interpolation method: ','s');

 disp(' ')
 descript = input('Type descriptor: ','s');

 disp(' ')
 fname = input('Type name of output file (without ext): ','s');

 Ofile = ['BCs/' fname '.mat'];

 switch Copt
 case {1,2}
   varout = 'descript t t_units X Y BCtype BC Imethod';
 case {3,4}
   varout = 'descript t t_units Y Z BCtype BC Imethod';
 case {5,6}
   varout = 'descript t t_units X Z BCtype BC Imethod';
 end

 eval(['save ' Ofile ' ' varout]);
 