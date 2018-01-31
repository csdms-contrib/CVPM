 function QS = QSsub_Z(S_opt,S0,hs,zf)

% Initializes the integral of the heat source over the control volumes
% (vertical dimension).

% Currently available heat sources:

%   S_opt = 'zero'           no heating
%         = 'linear'         source decreases linearly with depth
%         = 'exponential'    source decreases exponentially with depth

% To allow for a variable lithologic profile, different layers are allowed to 
% have different S0 and hs values.

% Notation:

%   S_opt = source function option                  (character string)
%   S0    = source function at Z = 0                (M+1)
%   hs    = source function length scale            (M+1)
%   zf    = CV interfaces                           (M+1)
%   QS    = source function integral                (M+1)
% ______________________________________________

% Written by:

%   Gary Clow
%   Institute of Arctic and Alpine Research
%   University of Colorado
%   Boulder, Colorado USA
%   Email: gary.clow@colorado.edu
% ______________________________________________

 M = length(zf) - 1;

% pre-allocate array

 QS = NaN*ones(size(zf));

% integrate S(z) over the control volumes

 switch S_opt
 case 'zero'
   QS = zeros(size(zf));

 case 'linear'
   for k=2:M
     QS(k) = S0(k) * ((zf(k+1) - zf(k)) - (zf(k+1)^2 - zf(k)^2)/(2*hs(k)));
   end

 case 'exponential'
   for k=2:M
     QS(k) = S0(k)*hs(k) * (expm1(-zf(k)/hs(k)) - expm1(-zf(k+1)/hs(k)));
   end
 end
