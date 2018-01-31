 function Kp = K_pore_3phase(Ki,Ku,Ka,psi_i,psi_u,psi_a)

% Finds the thermal conductivity of pores within rocks and soils using the
% Brailsford and Major 3-phase mixing model.  The pores are assumed to 
% consist of a mixture of ice, unfrozen water, and air.

% Notation:

%   Ki    = conductivity of ice                         (vector)
%   Ku    = conductivity of unfrozen water              (vector)
%   Ka    = conductivity of air                         (vector)
%   psi_i = relative volume fraction of ice             (vector)
%   psi_u = relative volume fraction of unfrozen water	(vector)
%   psi_a = relative volume fraction of air             (vector)
%   Kp    = conductivity of the pore space              (vector)
% ______________________________________________

% Written by:

%   Gary Clow
%   Institute of Arctic and Alpine Research
%   University of Colorado
%   Boulder, Colorado USA
%   Email: gary.clow@colorado.edu
% ______________________________________________

% parameters

 alpha = 0.75;
 f     = (1 - alpha) /2;
 g     = 2*alpha - 1;

% find conductivities assuming each phase is continuous

 Kic = K_BM3(Ki,Ku,Ka,psi_i,psi_u,psi_a);   % assume ice is continuous phase
 Kuc = K_BM3(Ku,Ki,Ka,psi_u,psi_i,psi_a);   % assume water is continuous phase
 Kac = K_BM3(Ka,Ki,Ku,psi_a,psi_i,psi_u);   % assume air is continuous phase

% find averaging weights

 wi = zeros(size(psi_i));
 wu = zeros(size(psi_i));
 wa = zeros(size(psi_i));
 N  = length(psi_i);

 for i=1:N

   if psi_i(i) >= alpha         % ice dominates
     wi(i) = 1;
     wu(i) = 0;
     wa(i) = 0;
   elseif psi_u(i) >= alpha     % unfrozen water dominates
     wu(i) = 1;
     wi(i) = 0;
     wa(i) = 0;
   elseif psi_a(i) >= alpha     % air dominates
     wa(i) = 1;
     wi(i) = 0;
     wu(i) = 0;
   else                         % no-man's land
     h      = (alpha - psi_i(i)) /2;    % no-man's wi(i)
     psi_up = psi_u(i) - h;
     psi_ap = psi_a(i) - h;
     if psi_ap <= f
       wi(i) = 1 - h / (g + psi_ap);
     else
       wi(i) = 1 - h / (g + psi_up);
     end

     h      = (alpha - psi_u(i)) /2;    % no-man's wu(i)
     psi_ip = psi_i(i) - h;
     psi_ap = psi_a(i) - h;
     if psi_ip <= f
       wu(i) = 1 - h / (g + psi_ip);
     else
       wu(i) = 1 - h / (g + psi_ap);
     end

     h      = (alpha - psi_a(i)) /2;    % no-man's wa(i)
     psi_ip = psi_i(i) - h;
     psi_up = psi_u(i) - h;
     if psi_up <= f
       wa(i) = 1 - h / (g + psi_up);
     else
       wa(i) = 1 - h / (g + psi_ip);
     end
   end
 end

% renormalize so sum(weights) = 1

 wsum = wi + wu + wa;
 wi   = wi ./ wsum;
 wu   = wu ./ wsum;
 wa   = wa ./ wsum;

% find weighted average conductivity

 Kp = wi.*Kic + wu.*Kuc + wa.*Kac;
