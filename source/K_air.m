 function Ka = K_air(T,planet)

% Finds the thermal conductivity of air (earth, 0.1 Pa) 
% or CO2 gas (mars, 0.6 kPa).

% Uses an interpolation table based on correlating equations from
% Stephan and Laesecke (1985) and Vesovic et al (1990).

% Range of validity is  70 K < Tk < 1000 K for air.
% Range of validity is 200 K < Tk < 1000 K for CO2 gas.

% Notation:

%   T  = temperature (C)                 (vector)
%   Ka = thermal conductivity of air     (vector)
% ______________________________________________

% Written by:

%   Gary Clow
%   Institute of Arctic and Alpine Research
%   University of Colorado
%   Boulder, Colorado USA
%   Email: gary.clow@colorado.edu
% ______________________________________________

 Tk = T + 273.15;

% table of [Tk, K_air, K_co2] values

 A = [100 0.0089086  NaN;      110 0.0098416  NaN; ...
      120  0.010764  NaN;      130  0.011676  NaN; ...
      140  0.012579  NaN;      150  0.013471  NaN; ...
      160  0.014354  NaN;      170  0.015226  NaN; ...
      180  0.016088  NaN;      190  0.016939  NaN; ...
      200  0.017779 0.0095617; 210  0.018609  0.010198; ...
      220  0.019428  0.010838; 230  0.020237  0.011496; ...
      240  0.021035  0.012178; 250  0.021822  0.012885; ...
      260    0.0226  0.013616; 270  0.023368  0.014368; ...
      280  0.024126   0.01514; 290  0.024874  0.015928; ...
      300  0.025614   0.01673; 310  0.026345  0.017544; ...
      320  0.027067  0.018367; 330  0.027782  0.019199; ...
      340  0.028488  0.020036; 350  0.029187  0.020879; ...
      360  0.029879  0.021725; 370  0.030564  0.022573; ...
      380  0.031242  0.023422; 390  0.031915  0.024272; ...
      400  0.032581  0.025122];

 Tk_air = A(:,1);

 switch planet
 case 'earth'
   Kair = A(:,2);
 case 'mars'
   Kair = A(:,3);
 end

 Ka = interp1(Tk_air,Kair,Tk);

 switch planet
 case 'earth'
   if any(Tk < 70)
     disp(' ')
     disp('K_air: temperature is beyond range of validity.')
   end
 case 'mars'
   if any(Tk < 200)
     disp(' ')
     disp('K_air: temperature is beyond range of validity.')
   end
 end
