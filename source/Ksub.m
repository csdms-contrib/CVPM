 function K = Ksub(T,Mtyp,Km0,phi,phi_i,phi_u,planet)

% Finds the bulk thermal conductivity for each CV.  For convenience, K
% is also defined on the boundaries where it is set equal to that of the 
% adjacent CVs.

% Notation:

%   T     = temperature (C)                         (M+1)
%   Mtyp  = material type                           (M+1)
%   Km0   = matrix conductivity at 0 C              (M+1)
%   phi   = porosity                                (M+1)
%   phi_i = volume fraction of ice                  (M+1)
%   phi_u = volume fraction of unfrozen water       (M+1)
%   K     = bulk thermal conductivity               (M+1)
% ______________________________________________

% Written by:

%   Gary Clow
%   Institute of Arctic and Alpine Research
%   University of Colorado
%   Boulder, Colorado USA
%   Email: gary.clow@colorado.edu
% ______________________________________________

 M = length(T) - 1;

% pre-allocate arrays

 K = NaN*ones(size(T));

% find which material types are present

 [Mtypes,m] = which_Mtypes(Mtyp);

% find K for each material present

 for j=1:m
   L = Mtyp == Mtypes(j);

   switch fix(Mtypes(j))
   case 1                       % fixed properties (used for testing)

     K(L) = Km0(L);

   case 2                       % linerized ice (used for testing)

     K(L) = K_linIce(T(L),Km0(L));

   case 3                       % ice

     K(L) = K_ice(T(L));

   case {4,5,6,7,8,9,10,11,12,13,14,15} % rocks and soils

     K(L) = K_rock(T(L),Mtyp(L),Km0(L),phi(L),phi_i(L),phi_u(L),planet);

   case {20,21,22,23,24,25}     % organic-rich material

     K(L) = K_veg(T(L),Mtyp(L),Km0(L),phi(L),phi_i(L),phi_u(L),planet);

   case {30,31,32,33,34,35}     % borehole fluids

     Ftype = Mtypes(j) - 29;
     K(L)  = K_fluid(T(L),Ftype);

   case {40,41,42,43,44}        % metals

     K(L) = K_metal(T(L),Mtyp(L));
   end
 end
