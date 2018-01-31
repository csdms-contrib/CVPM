 function [phi_u,dphiudT] = phiu_tableI(T,Mtyp,Mw,phi_u_tab,T_tab)

% Interpolates the previously determined (phi_u_tab,T_tab) table to find
% the value of phi_u at temperature T for each CV.  The derivative
% dphi_u/dT is also found.

% Notation:

%   T         = temperature of each CV                          (M+1)
%   Mtyp      = material type                                   (M+1)
%   Mw        = mass of water (both phases) per unit volume     (M+1)
%   phi_u_tab = table: volume fraction of unfrozen water        (M+1,nT)
%   T_tab     = table: temperatures corresponding to phi_u_tab  (M+1,nT)
%   phi_u     = volume fraction of unfrozen water at T          (M+1)
%   dphiudT   = dphi_u/dT at T                                  (M+1)

% Notes:

%   (1) The (phi_u_tab,T_tab) table was constructed using the initial 
%       constraint phi_u <= initial(phi_w).  As time proceeds, max
%       phi_u values may be somewhat smaller.  This will be the case
%       when the profile warms and ice converts to unfrozen water.  Thus,
%       we need to impose a more stringent constraint on phi_u as the
%       temperature field changes.  This constraint is provided by the 
%       total mass of water per unit volume (Mw) which is presumed to 
%       remain constant.  Then, phi_u <= (Mw/rhou) at all times.
%   (2) If any(Mw) < -100, then it is assumed that this routine has been
%       called by CPS, in which case the constraint phi_u <= initial(phi_w)
%       is fine.  This constraint is provided by routine phiu_pore.
% ______________________________________________

% Written by:

%   Gary Clow
%   Institute of Arctic and Alpine Research
%   University of Colorado
%   Boulder, Colorado USA
%   Email: gary.clow@colorado.edu
% ______________________________________________

% parameters

 rhou = 1000;       % density of unfrozen water
 dT   = 1e-06;      % temperature step for finding dphi_u/dT

% pre-allocate arrays

 phi_u   = zeros(size(Mtyp));
 dphiudT = zeros(size(Mtyp));

% set flag indicating whether the routine has been called within CPS

 if any(Mw < -100)
   flag = 1;
 else
   flag = 0;
 end

% step through layers

 for k=1:length(Mtyp)

   switch fix(Mtyp(k))
   case 1                           % fixed properties (used for testing)
     phi_u(k)   = 0;
     dphiudT(k) = 0;

   case {2,3}                       % pure ice
     phi_u(k)   = 0;
     dphiudT(k) = 0;

   case {4,5,6,7,8,9,10,11,12,13,14,15,20,21,22,23,24,25}   % rocks and organic-rich materials

%   interpolate to find phi_u(T)

     phi_u(k) = interp1(T_tab(k,:),phi_u_tab(k,:),T(k));      % phi_u(T)
     phi_up   = interp1(T_tab(k,:),phi_u_tab(k,:),T(k)+dT);   % phi_u(T+dT)
     phi_um   = interp1(T_tab(k,:),phi_u_tab(k,:),T(k)-dT);	  % phi_u(T-dT)

%   apply maximum value constraint

     if ~flag
       maxphi_u = Mw(k) / rhou;

       if phi_u(k) > maxphi_u
         phi_u(k)  = maxphi_u ;
       end
       if phi_up > maxphi_u
         phi_up  = maxphi_u ;
       end
       if phi_um > maxphi_u
         phi_um  = maxphi_u ;
       end
     end

% find the derviative dphi_u/dT

     dphiudT(k) = (phi_up - phi_um) / (2*dT);
   end
 end
