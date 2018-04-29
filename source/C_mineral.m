 function rhocp = C_mineral(T,Mtyp,rhom,cpm0)

% Finds the volumetric heat capacity of mineral grains at temperature
% T, where cpm0 is the value at the reference temperature (20 C).
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

% Notation:

%   T     = temperature (C)                     (vector)
%   Mtyp  = material type                       (vector)
%   rhom  = mean density of mineral grains      (vector)
%   cpm0  = specific heat of minerals (20 C)    (vector)
%   rhocp = volumetric heat capacity            (vector)

% Available material types (Mtyp):

%   igneous/metamorphic rocks
%       4   quartz dominated
%       5   feldspar dominated
%       6   mica dominated
%       7   pyroxene and amphibole dominated
%       8   olivine dominated
%   sedimentary rocks
%       10  sandstones
%       11  mudrocks
%       12  carbonates
%       13  cherts
% ______________________________________________

 Tk = T + 273.15;       % kelvin temperature
 Tr = Tk / 298.15;      % reduced temperature

% > Define the shape function

 sf = ones(size(cpm0));

% quartz (sandstones, cherts)

 L = (Mtyp == 4) | (Mtyp == 10) | (Mtyp == 13);

 if any(L)
   a = [-0.057 1.398 -0.340];
   sf(L) = a(1) + a(2)*Tr + a(3)*Tr.^2;
 end

% feldspars

 L = Mtyp == 5;

 if any(L)
   a = [-0.030 1.459 -0.427];
   sf(L) = a(1) + a(2)*Tr + a(3)*Tr.^2;
 end

% micas

 L = Mtyp == 6;

 if any(L)
   a = [-0.116 1.589 -0.471];
   sf(L) = a(1) + a(2)*Tr + a(3)*Tr.^2;
 end

% pyroxenes and amphiboles

 L = Mtyp == 7;

 if any(L)
   a = [-0.201 1.746 -0.543];
   sf(L) = a(1) + a(2)*Tr + a(3)*Tr.^2;
 end

% olivine

 L = Mtyp == 8;

 if any(L)
   a = [-0.259 1.854 -0.593];
   sf(L) = a(1) + a(2)*Tr + a(3)*Tr.^2;
 end

% clay minerals (mud rocks)

 L = Mtyp == 11;

 if any(L)
   a = [-0.220 1.711 -0.491];
   sf(L) = a(1) + a(2)*Tr + a(3)*Tr.^2;
 end

% cabonate minerals (limestones, dolomites)

 L = Mtyp == 12;

 if any(L)
   a = [ 0.124 1.252 -0.376];
   sf(L) = a(1) + a(2)*Tr + a(3)*Tr.^2;
 end

% > Find the specific heat cp(T) = cp(0) * sf

 cpm = cpm0 .* sf;

 if any(Tk < 150) || any(Tk > 400)
   disp(' ')
   disp('C_mineral: temperature is beyond range of validity.')
 end

% > Find the volumetric heat capacity of the mineral grains

 rhocp = rhom .* cpm;
