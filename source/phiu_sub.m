 function phi_u = phiu_sub(phi,lambda,r,DeltaT,PC_opt)

% Finds the volume fraction of unfrozen water for matrix particles
% with an effective radius r based on the Cahn et al. (1992) model.

% Notation:

%   phi    = porosity                                   (vector)
%   lambda = interfacial melting parameter [K^1/3 m]    (vector)
%   r      = effective radius of matrix particles [m]   (vector)
%   DeltaT = (Tf - T)                                   (vector)
%   PC_opt = packing coefficient option                 (string)
%   phi_u  = volume fraction of unfrozen water          (vector)

% Notes:

%   PC_opt = 'sc' or 'fcc'
% ______________________________________________

% Written by:

%   Gary Clow
%   Institute of Arctic and Alpine Research
%   University of Colorado
%   Boulder, Colorado USA
%   Email: gary.clow@colorado.edu
% ______________________________________________

% parameters

 xi = 0.0259e-06;   % curvature coefficient [K m]

 switch PC_opt
 case 'sc'          % simple cubic packing (sc)
   a1 = 1.893;
   a2 = 3.367;
 case 'fcc'         % cubic close packing (fcc)
   a1 = 2.450;
   a2 = 8.572;
 otherwise
   disp(' ')
   disp('warning phiu_sub: unknown packing option.')
   pause
 end

% > Find phi_u

% reset DeltaT values corresponding to temperatures warmer than Tf to yield
% phi_u values equal to phi 

 L         = DeltaT <= 0;
 DeltaT(L) = 1e-06;

% term1: water film on surface of spheres and at grain boundaries

 t1    = lambda ./ (r .* DeltaT.^(1/3));
 phiu1 = a1 * t1;

% term2: water in crevices and grain-boundary edges

 t2    = xi ./ (r .* DeltaT);
 phiu2 = a2 * t2.^2;

% add terms

 phi_u = phiu1 + phiu2;

% phi_u cannot exceed phi

 L = phi_u > phi;        
 if length(phi) == 1
   phi_u(L) = phi;
 elseif length(phi) == length(phi_u)
   phi_u(L) = phi(L);
 else
   disp(' ')
   disp('phiu_sub: warning, unknown case.')
   pause
 end
