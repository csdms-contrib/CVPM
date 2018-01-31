 function T = initTz_numerSS(Ts,qb,Ke,dz,QS)

% Finds initial temperature field (vertical dimension) assuming steady-state
% conditions.

% Notation:

%   Ts   = surface temperature at the start_time        (scalar)
%   qb   = basal heat flux at the start_time            (scalar)
%   Ke   = effective conductivity                       (M+1)
%   dz   = distance between grid points                 (M+1)
%   QS   = source term integral                         (M+1)
%   T    = temperature at the CV grid points            (M+1)

% Notes: 
% (1) The heat flux computed at the upper boundary J(2) is correct.  
%     However, to find the temperature for the grid point in the upper CV, 
%     what we need is the flux intermediate between this grid point and the
%     upper boundary.  This will yield a more accurate value for temperature
%     T(2).
% (2) This code currently assumes upperBC_type = 'T' and lowerBC_type = 'q'.
% ______________________________________________

% Written by:

%   Gary Clow
%   Institute of Arctic and Alpine Research
%   University of Colorado
%   Boulder, Colorado USA
%   Email: gary.clow@colorado.edu
% ______________________________________________

 M    = length(Ke) - 1;
 Kedz = Ke ./ dz;

% pre-allocate arrays

 T = NaN*ones(size(Ke));
 J = zeros(size(Ke));

% find heat flux at each CV interface

 J(M+1) = qb;

 for k=M:-1:2
   J(k) = J(k+1) - QS(k);
 end

% estimate the heat flux between points U and P for the upper CV, store in J(2) 

 J(2) = (3*J(2) + J(3)) / 4;

% find temperature at CV grid points using composite media formula

 T(1) = Ts;

 for k=2:M+1
   T(k) = T(k-1) - J(k) / Kedz(k); 
 end
