 function K = K_veg(T,Mtyp,Km0,phi,phi_i,phi_u,planet)

% Finds the bulk thermal conductivity of organic-rich material at 
% temperature T.  The material is assumed to consist of a matrix of 
% organic and mineral particles while the pores are filled with liquid 
% water, ice, and air.

% Notation:

%   T     = temperature (C)                         (vector)
%   Mtyp  = material type                           (vector)
%   Km0   = mineral conductivity at 0 C             (vector)
%   phi   = porosity                                (vector)
%   phi_i = volume fraction of ice                  (vector)
%   phi_u = volume fraction of unfrozen water       (vector)
%   K     = bulk thermal conductivity               (vector)

% Note:
%   (1) This routine is identical to K_rock except the matrix has both
%       organic and mineral particles.
% ______________________________________________

% Written by:

%   Gary Clow
%   Institute of Arctic and Alpine Research
%   University of Colorado
%   Boulder, Colorado USA
%   Email: gary.clow@colorado.edu
% ______________________________________________

% pre-allocate arrays

 K = NaN*ones(size(phi));

% find the conductivity of the organic matrix

 Km = K_organic(T,Mtyp,Km0);

% find the bulk conductivity of porous material (use same code as for soils)

 K = K_rock_porous(T,Km,phi,phi_i,phi_u,planet);
