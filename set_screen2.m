 function [pos] = set_screen2(Ssize)

% Sets position and size of matlab plot window.
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

% set Ssize = 0 full screen
%           = 1 reserved
%           = 2 reserved
%           = 3 for manuscript (2-column format)
%           = 4 for manuscript (full page width)
%           = 5 for manuscript (MET time series)
% ______________________________________________

 switch Ssize
 case 0
   scn = get(0,'ScreenSize');
 case {1,2,3,4}
   scn = [1 1 1280 1024];   % use this as standard screen size
 end

 h   = scn(4);      % height
 w   = scn(3);      % width
 x0  = 0.25*scn(3); % x-coordinate lower left corner
 y0  = 0.25*scn(4); % y-coordinate lower left corner

 switch Ssize
 case 0
   rx = 1;                          % x-reduction factor
   ry = 1;                          % y-reduction factor
 case 3
   rx = 0.40;
   ry = 0.40;
   set(0,'DefaultAxesFontSize',16)	% axes labels
   set(0,'DefaultTextFontSize',16)	% text
 case 4
   rx = 0.7;
   ry = 0.7;
   set(0,'DefaultAxesFontSize',16)	% axes labels
   set(0,'DefaultTextFontSize',16)	% text
 case 5
   rx = 1.0;
   ry = 0.1;
   set(0,'DefaultAxesFontSize',16)	% axes labels
   set(0,'DefaultTextFontSize',16)	% text
 end
 pos = [x0 y0 (x0 + rx*0.50*w) (y0 + ry*0.50*h)];
