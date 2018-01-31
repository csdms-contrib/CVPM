 function K = K_BM3(K1,K2,K3,psi1,psi2,psi3)

% Finds the bulk thermal conductivity of a 3-component system using the 
% Brailsford and Major model that assumes phases 2 and 3 are randomly
% distributed within a continuous phase 1 medium.

% Notation:

%   K1   = conductivity of phase 1                 (vector)
%   K2   = conductivity of phase 2                 (vector)
%   K3   = conductivity of phase 3                 (vector)
%   psi1 = relative volume fraction of phase 1     (vector)
%   psi2 = relative volume fraction of phase 2     (vector)
%   psi3 = relative volume fraction of phase 3     (vector)
%   K    = bulk thermal conductivity               (vector)

% Notes:
%   (1)   1 = psi1 + psi2 + psi3      (this should be true)
% ______________________________________________

% Written by:

%   Gary Clow
%   Institute of Arctic and Alpine Research
%   University of Colorado
%   Boulder, Colorado USA
%   Email: gary.clow@colorado.edu
% ______________________________________________

 chi2 = (K1 ./ K2);
 chi3 = (K1 ./ K3);

 f2   = psi2 ./ (2*chi2 + 1);
 f3   = psi3 ./ (2*chi3 + 1);
 f4   = f2 + f3;
 f5   = f2 .* chi2 + f3 .* chi3;

 K    = K1 .* (psi1 + 3 * f4) ./ (psi1 + 3 * f5);
