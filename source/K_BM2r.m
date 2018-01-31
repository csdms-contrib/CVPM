 function K = K_BM2r(K1,K2,psi2)

% Finds the bulk thermal conductivity of a 2-component system using the 
% Brailsford and Major model that assumes both phases are randomly
% distributed (their eq. 13, corrected).

% Notation:

%   K1   = conductivity of phase 1                  (vector)
%   K2   = conductivity of phase 2                  (vector)
%   psi2 = relative volume fraction of phase 2      (vector)
%   K    = bulk thermal conductivity                (vector)

% Notes:
%   (1)     1 = psi1 + psi2
% ______________________________________________

% Written by:

%   Gary Clow
%   Institute of Arctic and Alpine Research
%   University of Colorado
%   Boulder, Colorado USA
%   Email: gary.clow@colorado.edu
% ______________________________________________

 chi = (K1 ./ K2);

 f1  = 2*chi - 1;
 f2  = 3*psi2 .* (chi - 1);
 f   = f1 - f2;

 K   = K1 .* (f + sqrt(f.^2 + 8*chi)) ./ (4*chi);
