 function [Jtilde,Jnet] = Jsub_R(BCtypes,T,rf,Ke,dr,qo)

% Finds the integrated heat fluxes across CV interfaces and the net
% heat flux out of each control volume.  J and Jnet calculated by
% this function are not used to find the temperature field.  Instead,
% they are provide to assist with CVPM diagnostics.

% Notation:

%   T      = temperature                              (N+1)
%   rf     = interface locations                      (N+1)
%   Ke     = effective conduct. across R-interfaces   (N+1)
%   qo     = heat flux through outer interface        (scalar)
%   Jtilde = integrated heat flux at CV interfaces    (N+1)
%   Jnet   = net heat-flux out of each CV             (N+1)

% Notes: 
%   (1) qb is not used if lowerBC_type = 'T'.
%   (2) For the inner CV, the heat flux J(2) initially calculated 
%       from the temperature difference T(2)-T(1) pertains to a point
%       midway between grid points W and P.  This value and J(3) are
%       are then used to find the flux at the inner boundary.  This
%       is then stored in J(2).
% ______________________________________________

% Written by:

%   Gary Clow
%   Institute of Arctic and Alpine Research
%   University of Colorado
%   Boulder, Colorado USA
%   Email: gary.clow@colorado.edu
% ______________________________________________

 N = length(T) - 1;
 innerBC_type = BCtypes{1};
 outerBC_type = BCtypes{2};

% pre-allocate arrays

 Jtilde = NaN*ones(size(T));
 Jnet   = NaN*ones(size(T));

% find the heat flux across each CV interface

 rKedr = rf .* Ke ./ dr;

 for k=2:N
   Jtilde(k) = -rKedr(k) * (T(k) - T(k-1));
 end

% find an improved value for the flux on the inner boundary

 switch innerBC_type
 case 'T'
   Jtilde(2) = (4*Jtilde(2) - Jtilde(3)) / 3;
 end

 switch outerBC_type
 case 'q'
   Jtilde(N+1) = 2*pi*rf(N+1)*qo;
 end

% find net heat-flux out of each CV

 for k=2:N
   Jnet(k) = Jtilde(k+1) - Jtilde(k);
 end
