 function [td,Ta] = SB_drill_model(dR,Rdo,Rh,H,rho,mf,Te,T0,gamma,K,kappa,tdays,Z)

% Finds temperatures inside a drill pipe and surrounding annulus using the  
% Szarka-Bobok model,

% Szarka, Z., and Bobok, E. (2012): Determination of the temperature distribution
%   in the circulating drilling fluid, Geosci. Eng., 1, 37-47.
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

% This version is designed to simulate the advancing drill.

% Notation:

%   dR    = drill pipe wall thickness [m]
%   Rdo   = drill pipe outer radius [m]
%   Rh    = open hole radius [m]
%   H     = total hole depth [m]
%   rho   = fluid density [kg/m^3]
%   mf    = fluid mass flux [kg/s]
%   Te    = fluid inlet temperature [C]
%   T0    = mean long-term surface temperature [C]
%   gamma = mean geothermal gradient [K/m]
%   K     = mean thermal conductivity of rock
%   kappa = mean thermal diffusivity of rock
%   tdays = time to drill hole [days]
%   Z     = depth grid [m]
%   td    = time drill bit reaches each depth in vector Z [days]
%   Ta    = annulus temperatures [C], Ta(:,i) is the T-profile at time td(i)

% Options:

 dflag = 0;   % set dflag = 1 to show diagnostics
 mflag = 0;	  % set mflag = 1 to simulate presence of mudcake at borehole wall
% ______________________________________________

% > Parameters

% pipe properties

 Kd = 45;                       % thermal conductivity of drill pipe [W/(m K)]

% fluid properties

 c       = 0.9*4180;            % heat capacity of fluid [J/(kg K)]
 Kf      = 0.77;                % thermal conductivity [W/(m K)];
 kappaf  = 1.5e-07;             % thermal diffusivity [m^2/s]
 nu      = 1.1e-05;             % kinematic viscosity [m^2/s]

% time used in transient heat-conduction function

 t_hcf = 86400 * 0.01;          % i.e., 0.01 day

% > Begin calculations

% undisturbed formation temperatures

 F = T0 + gamma*Z;

% mean drilling velocity

 Vd = H / tdays;    % [m/day]

% find time drill bit reaches each depth Z

 td = Z' / Vd;  % [days}

% find inner radius of drill pipe

 Rdi = Rdo - dR;

% set inner radius of casing

 Rci = Rh;

% find drill-fluid velocities

 Vf = mf / rho;                 % fluid volume flux [m^3/s]
 Ad = pi * Rdi^2;               % cross sectional area inside drill pipe [m^2]
 vd = Vf / Ad;                  % fluid velocity inside drill pipe [m/s]
 Aa = pi * (Rci^2 - Rdo^2);     % cross sectional area of annulus
 va = Vf / Aa;                  % fluid velocity inside annulus

% Prandtl number

 Pr = nu / kappaf;

% Reynolds number

 Re_d = vd * (2*Rdi) / nu;                  % in drill pipe
 Re_a = va * (2*sqrt(Rci^2 - Rdo^2)) / nu;  % in annulus

% Nusselt number

 Nu_d = 0.027 * Re_d^(4/5) * Pr^(1/3);      % in drill pipe
 Nu_a = 0.027 * Re_a^(4/5) * Pr^(1/3);      % in annulus

% heat transfer coeffients

 h1 = Nu_d * Kf / (2*Rdi);
 h2 = Nu_a * Kf / (2*Rdo);

% drillpipe/borehole parameters

 x    = 1/h1 + dR/Kd;
 Udi  = 1 / x;
 x    = 1/h2 + dR/Kd;
 Udo  = 1 / x;
 x    = 1/h2 + dR/K;
 Uci1 = 1 / x;          % with mudcake
 Uci2 = h2;             % without mudcake
 if mflag
   Uci = Uci1;
 else
   Uci = Uci2;
 end

 if dflag
   disp(' ')
   disp(['Re_d = ' num2str(Re_d)])
   disp(['Re_a = ' num2str(Re_a)])
   disp(['Nu_d = ' num2str(Nu_d)])
   disp(['Nu_a = ' num2str(Nu_a)])
   disp(['h1   = ' num2str(h1)])
   disp(['h2   = ' num2str(h2)])
   disp(['Udi  = ' num2str(Udi)])
   disp(['Udo  = ' num2str(Udo)])
   disp(['Uci  = ' num2str(Uci1) '  (with mudcake on outer wall)'])
   disp(['Uci  = ' num2str(Uci2) '  (without mudcake)'])
 end

% transient heat-conduction function

 x = kappa * t_hcf / Rdo^2;
 y = log(sqrt(4*x));
 p = [0.08426 0.51066 0.62076];
 f = polyval(p,y);

% borehole performance state coef [m]

 fac1 = Rci*Uci*f;
 fac2 = 2*pi*Rci*Uci*K;
 A    = mf*c*(K + fac1)/fac2;

% drillpipe performance state coef [m]

 B = (mf*c) / (2*pi*Rdi*Udi);

% roots of characteristic equation

 fac3    = sqrt(1 + 4*A/B);
 lambda1 = (1 + fac3)/(2*A);
 lambda2 = (1 - fac3)/(2*A);

% > Time loop

 nZ = length(Z);
 nt = length(td);
 T  = NaN*ones(nZ,nt);

 for i=1:nt

   H = Z(i);    % depth of hole at time td(i)

% find various constants

   D  =  lambda2*exp(lambda2*H) - lambda1*exp(lambda1*H);
   D1 =  lambda2*(Te - T0 + B*gamma)*exp(lambda2*H) + gamma;
   D2 = -lambda1*(Te - T0 + B*gamma)*exp(lambda1*H) - gamma;
   C1 =  D1/D;
   C2 =  D2/D;

% find temperatures in the annulus

   fac5      = NaN*ones(size(Z));
   fac5(1:i) = C1*(1+B*lambda1)*exp(lambda1*Z(1:i)) + C2*(1+B*lambda2)*exp(lambda2*Z(1:i));
   Ta        = fac5 + F;
   T(:,i)    = Ta; 
 end

% Correct a small computational error in the upper couple meters during the 
% first couple of time steps.

 for k=1:4
   for i=1:5
     if T(k,i) > Te
       T(k,i) = Te;     % T(k,i) can't be greater than Te here
     end
   end
 end

 Ta = T;
