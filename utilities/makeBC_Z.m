% makeBC_Z.m

% Make a boundary-condition file for the 1-D vertical CVPM case.
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

% Notes:

%   (1) At present, this code is only setup to create an upper BC using  
%       North Slope air temperature records as a proxy for the upper BC
%       (surface temperature).
%   (2) For the air temperature, the Bieniek air-temp record is welded to
%       the local DOI/GTN-P air-temp record as year 1999.
%   (3) The composite air-temp record is then smoothed by passing it through 
%       a N-yr running average.  This smoothed record is then offset to 
%       produce a proxy for the surface ground-temperature record.
% __________________________________________________________________

% > Set environment

 close all
 clear all
 format shortg
 colordef white
 pos = set_screen(0);

% set text and symbol sizes

 set(0,'DefaultAxesFontSize',20)

% options:

 Nyr   = 9;     % number of years in moving average
 tflag = 0;     % tflag=1 resets first year to zero for output time-series

% set working directory

 wdir = '~/thermal/numer/CVPM';
 cd(wdir)

% set directory containing the MET series

 mdir = '~/thermal/numer/CVPM/AK_projects/MET';
% ________________________________________________________________

 disp(' ')
 disp(['*** This code currently set to use ' num2str(Nyr) ' year averaging. ***'])

 disp(' ')
 well_code = input('Type well_code: ','s');
 disp(' ')
 Ts0_lo = input('Type cold end of long-term Ts from extrapolated borehole temps: ');
 Ts0_hi = input('Type warm end of long-term Ts from extrapolated borehole temps: ');

 Ts0  = (Ts0_lo:0.1:Ts0_hi);
 nTs0 = length(Ts0);

 disp(' ')
 disp('Available Air Temperature Series:')
 disp('[1] Chukchi coast    (TLK) + Bieniek')
 disp('[2] Beaufort coast   (DRP) + Bieniek')
 disp('[3] NE coastal plain (FCK) + Bieniek')
 disp('[4] Arctic foothills (UMI) + Bieniek')
 Aopt = input('Select series (1,2,...): ');

% > Load DOI/GTN-P temperature series

 switch Aopt
 case 1
   Ifile = 'region_Chukchi_series';
   s1    = ' DOI/GTN-P Chukchi series';
 case 2
   Ifile = 'region_Beaufort_series';
   s1    = ' DOI/GTN-P Beaufort series';
 case 3
   Ifile = 'region_Northeast_series';
   s1    = ' DOI/GTN-P NE-ACP series';
 case 4
   Ifile = 'region_Foothills_series';
   s1    = ' DOI/GTN-P Foothills series';
 end

 load([mdir '/' Ifile])
 tyr       = ann_seriesC{1};
 Tair_ann  = ann_seriesC{2};
 Tg10_ann  = ann_seriesC{3};
 Tg120_ann = ann_seriesC{4};

% find regional average temperatures and offsets over (2000-2008)

 Tair  = Tair_ann;
 Tg10  = Tg10_ann;
 Tg120 = Tg120_ann;

 L     = tyr >= 2000 & tyr <= 2008;
 Tair  = Tair(L);
 Tg10  = Tg10(L);
 Tg120 = Tg120(L);

 Tair_ave  = ave_static_nan(Tair_ann,0.1);
 Tg10_ave  = ave_static_nan(Tg10_ann,0.1);
 Tg120_ave = ave_static_nan(Tg120_ann,0.1);

 Tg10_offset  = Tg10_ave  - Tair_ave;
 Tg120_offset = Tg120_ave - Tair_ave;

 Tair_ave
 Tg10_ave
 Tg120_ave
 Tg10_offset
 Tg120_offset

 clear Tair Tg10 Tg120 Tair_ave Tg10_ave Tg120_ave

% > Load Bieniek air-temp anomaly record

 Bfile = '2014Bieniek_ann_anom.mat';

 load([mdir '/' Bfile])
 Byr     = Acell{1};
 B_Tanom = Acell{3};

% strip off bieniek 2012 which is bogus (doesn't match DRP)

 L       = Byr <= 2011;
 Byr     = Byr(L);
 B_Tanom = B_Tanom(L);
 nB      = length(Byr);

% find offset between the Bieniek record and our air-temp record

 L1 = tyr >= 1999 & tyr <= 2011;
 x  = Tair_ann(L1);
 L2 = Byr >= 1999 & Byr <= 2011;
 y  = B_Tanom(L2);

 B_offset    = y - x;
 Boffset_ave = mean(B_offset);
 Boffset_std = std(B_offset);

 Boffset_ave
 Boffset_std

% find Bieniek air temps adjusted to this particular site

 B_Tair = B_Tanom - Boffset_ave;

% weld the Bieniek air temps and our air temps together at year 1999

 L     = Byr < 1999;
 year  = Byr(L);
 Tair  = B_Tair(L);
 L     = tyr >= 1999 & tyr <= 2015;
 year  = [year; tyr(L)];
 Tair  = [Tair; Tair_ann(L)];

% create smoothed versions of Tair

 TairN  = ave_moving(Tair,Nyr);     % N-yr moving average
 TairN3 = ave_moving(TairN,3);      % remove remaining spikes with a 3-pt moving ave

% > Create a ground-temperature series from the smoothed air-temperature series

 disp(' ')
 disp('The BC will be based on the local Tair series, offset by an appropriate amount.')
 disp(' ')
 Tgoffset_lo = input('Type lower end of Tgoffset range (eg., 2.6): ');
 Tgoffset_hi = input('Type upper end of Tgoffset range (eg., 3.1): ');

 Tgoffset = (Tgoffset_lo:0.1:Tgoffset_hi);
 noffset  = length(Tgoffset);

 for j=1:nTs0
   for k=1:noffset

     Ts = TairN3 + Tgoffset(k);    % Tg estimated from our Tair

     if j == 1
       h1 = figure('position',pos);
       ax(1) = subplot(2,1,1);
       plot(year,Tair,'linestyle','none','color','b','marker','o','markeredgecolor','b','markerfacecolor','y')
       hold on
       plot(year,TairN3,'color','b')
       grid on
       ylabel('$T$ ($^\circ$C)','interpreter','latex')
       title('Air Temperature','interpreter','latex')

       ax(2) = subplot(2,1,2);
       plot(year,Ts,'color','m','linewidth',2)
       hold on
       plot(tyr,Tg120_ann,'color','r','marker','s','markeredgecolor','r','markerfacecolor','y','linewidth',1)
       grid on
       linkaxes(ax,'x')
       xlabel('Year','interpreter','latex')
       ylabel('$T$ ($^\circ$C)','interpreter','latex')
       title(['Proxy Ground Temperature (120 cm), Tgoffset = ' num2str(Tgoffset(k))],'interpreter','latex')
       pause
     end

%   tack the long-term value on the front end

     t  = [1880; year];
     BC = [Ts0(j); Ts];      

%   reset first time to zero if desired

     if tflag
       t = t - t(1);
     end

     if BC(1) <= BC(2)      % if T(1880) <= T(1920), then show plot and output file

       h2 = figure('position',pos);
       plot(t,BC,'color','r','linewidth',2)
       grid on
       xlabel('Year','interpreter','latex')
       ylabel('$T$ ($^\circ$C)','interpreter','latex')
       title(['Boundary Condition, Ts0 = ' num2str(Ts0(j)) ', Tgoffset = ' num2str(Tgoffset(k))])
       if j == 1
         pause
       else
         pause(2)
       end
       close(h2)

%     setup information needed by CVPM

       descript = [s1 ', Tair offset by ' num2str(Tgoffset(k))];
       t_units  = 'years';
       BCtype   = 'T';
       Imethod  = 'linear';
       varout   = 'descript t t_units BCtype BC Imethod';

       s2 = num2str((round(10*abs(Ts0(j)))/10));
       ip = find(s2 == '.');
       if ~isempty(ip)
         s2(ip) = 'p';
       else
         s2 = [s2 'p0'];
       end

       s3 = num2str((round(1000*abs(Tgoffset(k)))/1000));
       ip = find(s3 == '.');
       s3(ip) = 'p';
       if ~isempty(ip)
         s3(ip) = 'p';
       else
         s3 = [s3 'p0'];
       end

       fname = ['Ts_' well_code '_m' s2 '_air' s3 '.mat'];
       Ofile = ['BCs/' fname];

       eval(['save ' Ofile ' ' varout]);
     end
   end
 end
