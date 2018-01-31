 function Km = K_mineral(T,Mtyp,Km0)

% Finds the thermal conductivity of mineral grains at temperature T.
% This routine takes into account the difference between 
% igneous/sedimentary and sedimentary mineral assemblages.

% Notation:

%   T    = temperature (C)                              (vector)
%   Mtyp = material type                                (vector)
%   Km0  = mineral conductivity at 0 C                  (vector)
%   Km   = mineral conductivity at temperature T        (vector)
% ______________________________________________

% Written by:

%   Gary Clow
%   Institute of Arctic and Alpine Research
%   University of Colorado
%   Boulder, Colorado USA
%   Email: gary.clow@colorado.edu
% ______________________________________________

% pre-allocate the array

 Km = NaN*ones(size(Km0));

% igneous/metamorphic minerals

 L = (Mtyp >= 4) & (Mtyp <= 9);

 if any(L)
   a = [0.99 0.0030 0.0042];
   Km(L) = Km0(L) ./ (a(1) + T(L) .* (a(2) - a(3)./Km0(L)));
 end

% sedimentary minerals

 L = (Mtyp >= 10) & (Mtyp <= 15);

 if any(L)
   a = [0.99 0.0034 0.0039];
   Km(L) = Km0(L) ./ (a(1) + T(L) .* (a(2) - a(3)./Km0(L)));
 end
