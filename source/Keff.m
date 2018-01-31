 function Ke = Keff(K,varep)

% Finds the effective conductivity at CV interfaces.

% Notation:

%   K     = conductivity at the CV grid pts             (M+1)
%   varep = fractional distance                         (M+1)
%   Ke    = effective conductivity at CV interfaces     (M+1)
% ______________________________________________

% Written by:

%   Gary Clow
%   Institute of Arctic and Alpine Research
%   University of Colorado
%   Boulder, Colorado USA
%   Email: gary.clow@colorado.edu
% ______________________________________________

 M = length(K) - 1;

% pre-allocate array

 Ke = NaN*ones(size(varep));

% find effective conductivity

 for k=3:M
   Ke(k) = 1 / ( (1-varep(k))/K(k-1) + varep(k)/K(k) );
 end

 Ke(1:2) = K(2);    % effective Ke on upper boundary
 Ke(M+1) = K(M);    %     "     "     lower boundary
