 function Ku = K_water(T)

% Finds the thermal conductivity of liquid water.
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

% Uses an interpolation table generated from Huber et al's (2012)
% correlating equation.
% Range of validity is 250 K < T < 383.15 K.
% For T < 250 K, the value at 250 K is assumed while for T > 383.15 K, the 
% value at 383.15 K is assumed.

% Notation:

%   T  = temperature (C)                            (vector)
%   Ku = thermal conductivity of unfrozen water     (vector)
% ______________________________________________

 Tk = T + 273.15;

% table of K_water(Tk) values

 A = [250 0.47247; 252 0.48236; 254 0.49155; 256 0.50011; ...
      258 0.50812; 260 0.51563; 262 0.52269; 264 0.52934; ...
      266 0.53563; 268 0.54159; 270 0.54726; 272 0.55266; ...
      274 0.55781; 276 0.56274; 278 0.56746; 280 0.57199; ...
      282 0.57635; 284 0.58055; 286  0.5846; 288 0.58852; ...
      290  0.5923; 292 0.59595; 294  0.5995; 296 0.60293; ...
      298 0.60626; 300 0.60949; 302 0.61262; 304 0.61566; ...
      306 0.61861; 308 0.62148; 310 0.62426; 312 0.62696; ...
      314 0.62959; 316 0.63213; 318  0.6346; 320   0.637; ...
      322 0.63932; 324 0.64158; 326 0.64376; 328 0.64588; ...
      330 0.64793; 332 0.64991; 334 0.65182; 336 0.65367; ...
      338 0.65546; 340 0.65718; 342 0.65884; 344 0.66044; ...
      346 0.66198; 348 0.66345; 350 0.66487; 352 0.66623; ...
      354 0.66754; 356 0.66878; 358 0.66997; 360  0.6711; ...
      362 0.67218; 364 0.67321; 366 0.67418; 368  0.6751; ...
      370 0.67597; 372 0.67678; 374 0.67755; 376 0.67827; ...
      378 0.67894; 380 0.67956; 382 0.68013; 384 0.68066];

 Tk_water = A(:,1);
 Ku_water = A(:,2);

% interpolate to find K_water

 Ku = interp1(Tk_water,Ku_water,Tk);

% deal with temperatures beyond range of validity

 if any(Tk < 250)
   L     = Tk < 250;
   Ku(L) = interp1(Tk_water,Ku_water,250);
 end

 if any(Tk > 383.15)
   L     = Tk > 383.15;
   Ku(L) = interp1(Tk_water,Ku_water,383.15);
 end
