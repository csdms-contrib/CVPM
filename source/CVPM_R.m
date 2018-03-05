 function [] = CVPM_R(experim)

% Computational engine for the CVPM 1-D radial case (R).
% ______________________________________________

% Written by:

%   Gary Clow
%   Institute of Arctic and Alpine Research
%   University of Colorado
%   Boulder, Colorado USA
%   Email: gary.clow@colorado.edu
% ______________________________________________

 disp(' ')
 disp('* Running CVPM *')
 disp(' ')

% import initial fields created by CPS

 Ifile = ['CPSout/' experim '_cps.mat'];

 load(Ifile)
 f = IEfac;

 disp(['dt = ' num2str(tstep) ' ' t_units ' (computational time step)'])
 disp(' ')

% unpack time limits

 tmin = t_limits(1);
 tmax = t_limits(2);

% get boundary condition arrays

 innerBC_type = BCtypes{1};
 outerBC_type = BCtypes{2};

 innerBC_file = ['BCs/' BCfiles{1}];
 outerBC_file = ['BCs/' BCfiles{2}];

 switch innerBC_type
 case 'T'
   [des,t_TaA,tunits_Ta,TaA,Tamethod] = inputBC(innerBC_file,CS);
   t_TaA = convert_tunits(innerBC_file,t_TaA,tunits_Ta,t_units);
   Ta    = interp1(t_TaA,TaA,tmin,Tamethod);
   qa    = 0;
 case 'q'
   [des,t_qaA,tunits_qa,gaA,qamethod] = inputBC(innerBC_file,CS);
   t_qaA = convert_tunits(innerBC_file,t_qaA,tunits_qa,t_units);
   qa    = interp1(t_qaA,qaA,tmin,qamethod);
   Ta    = 0;
 case 'none'
   Ta    = T(1);
   qa    = 0;
 end

 switch outerBC_type
 case 'T'
   [des,t_ToA,tunits_To,ToA,Tomethod] = inputBC(outerBC_file,CS);
   t_ToA = convert_tunits(outerBC_file,t_ToA,tunits_To,t_units);
   To    = interp1(t_ToA,ToA,tmin,Tomethod);
   qo    = 0;
 case 'q'
   [des,t_qoA,tunits_qo,qoA,qomethod] = inputBC(outerBC_file,CS);
   t_qoA = convert_tunits(outerBC_file,t_qoA,tunits_qo,t_units);
   qo    = interp1(t_qoA,qoA,tmin,qomethod);
   To    = 0;
 end

% setup computational time grid

 secs_tunit = dtunit2secs(t_units);     % seconds per t_unit

 dt    = secs_tunit * tstep;            % calculation interval (secs)
 Nt    = round((tmax-tmin)/tstep);      % number of time steps
 tgrid = tmin + tstep * (1:Nt);         % computational time grid (in t_units)

% setup output time grid

 if ostep < 0           % variable time step
   Ostep = abs(ostep);
   ogrid = [tmin (Ostep:Ostep:10*Ostep)];
   Ostep = 10*Ostep;
   ogrid = [ogrid (2*Ostep:Ostep:10*Ostep)];
   Ostep = 10*Ostep;
   ogrid = [ogrid (2*Ostep:Ostep:10*Ostep)];
   Ostep = 10*Ostep;
   ogrid = [ogrid (2*Ostep:Ostep:10*Ostep)];
   Ostep = 10*Ostep;
   ogrid = [ogrid (2*Ostep:Ostep:100*Ostep)];
   L     = ogrid <= tmax;
   ogrid = ogrid(L);
 else                   % uniform time step
   ogrid = (tmin:ostep:tmax);
 end
 Nout = length(ogrid);   % number of output times

% pre-allocate storage arrays

 tarr     = NaN*ones(  1,Nout);
 Taarr    = NaN*ones(  1,Nout);
 qaarr    = NaN*ones(  1,Nout);
 Toarr    = NaN*ones(  1,Nout);
 qoarr    = NaN*ones(  1,Nout);

 Tarr     = NaN*ones(N+1,Nout);
 phi_iarr = NaN*ones(N+1,Nout);
 phi_uarr = NaN*ones(N+1,Nout);
 Karr     = NaN*ones(N+1,Nout);
 Carr     = NaN*ones(N+1,Nout);
 Jarr     = NaN*ones(N+1,Nout);
 Jnetarr  = NaN*ones(N+1,Nout);
 QSarr    = NaN*ones(N+1,Nout);

% initial values for storage arrays

 Ke            = Keff(K,varepR);
 [Jtilde,Jnet] = Jsub_R(BCtypes,T,rf,Ke,dr,qo);

 istor             = 1;
 tarr(      istor) = tmin;
 Taarr(     istor) = Ta;
 qaarr(     istor) = qa;
 Toarr(     istor) = To;
 qoarr(     istor) = qo;
 Tarr(    :,istor) = T';
 phi_iarr(:,istor) = phi_i';
 phi_uarr(:,istor) = phi_u';
 Karr(    :,istor) = K';
 Carr(    :,istor) = C';
 Jarr(    :,istor) = Jtilde';
 Jnetarr( :,istor) = Jnet';

% find source-function integrated over the CVs (assume it's zero for the 1D radial case)

   QS = zeros(size(rf));

% ------ Begin Time Loop -------->

 Tn  = T;
 qan = qa;
 qon = qo;

 for i=1:Nt

   t = secs_tunit * tgrid(i);           % time (secs)

% find BCs

   switch innerBC_type
   case 'T'
     Ta = interp1(t_TaA,TaA,tgrid(i),Tamethod);
     qa = 0;
   case 'q'
     Ta = 0;
     qa = interp1(t_qaA,qaA,tgrid(i),qamethod);
   case 'none'
     Ta = Tn(1);
     qa = 0;
   end

   switch outerBC_type
   case 'T'
     To = interp1(t_ToA,ToA,tgrid(i),Tomethod);
     qo = 0;
   case 'q';
     To = 0;
     qo = interp1(t_qoA,qoA,tgrid(i),qomethod);
   end

   Tbc = [Ta  To];
   q   = [qa  qo];
   qn  = [qan qon];

% find the volume fractions of ice and unfrozen water

   [phi_u,dphiudT] = phiu_tableI(T,Mtyp,Mw,phi_u_tab,T_tab);
   phi_i           = phiiSub(Mtyp,Mw,phi_u);

% find the thermal conductivity at the CV grid points

   K = Ksub(T,Mtyp,Km0,phi,phi_i,phi_u,planet);

% find effective thermal conductivity at CV interfaces

   Ke = Keff(K,varepR);

% find the volumetric heat capacity

   [rhocp,C] = Csub(T,Mtyp,rhom,cpm0,phi,phi_i,phi_u,dphiudT);

% find discretization coefficients

   [aW,aWp,aE,aEp,aP,aPp,b] = DEcoefs_R(dt,f,rf,Ar,VP,C,Ke,QS,BCtypes,q,qn);

% find new temperature field

   T = TDMA_R(aW,aP,aE,aWp,aPp,aEp,b,Tn,BCtypes,Tbc);

% estimate temperature on outer boundary if outerBC_type = 'q'

   switch outerBC_type
   case 'q'
     T(N+1) = T(N) - qo * dr(N+1) / Ke(N+1);
   end

% diagnostics: find the integrated heat-flux across the interfaces

   [Jtilde,Jnet] = Jsub_R(BCtypes,T,rf,Ke,dr,qo);

% check numerical stability

   VPC      = VP .* C;
   AKeR     = Ar .* Ke;
   fac      = NaN*ones(1,N+1);
   fac(1:N) = AKeR(1:N) + AKeR(2:N+1);

   sfac     = VPC ./ ((1-f) * fac);

   switch innerBC_type
   case {'T','q'}
     jmin = 3;
   case 'none'
     jmin = 2;
   end
   minSfac = min(sfac(jmin:N));

   if minSfac < dt
     disp(' ')
     disp(['Warning: time step is too big for numerical stability, t = ' num2str(tgrid(i)) ' (' t_units ')'])
     secs_tunit = dtunit2secs(t_units);     % seconds per t_unit
     minSfac    = minSfac / secs_tunit;     % convert to t_units
     disp(' ')
     disp(['dt should be less than ' num2str(minSfac) ' (' t_units ') to guarantee numerical stability.'])
   end

% store results for output

   tfac = abs(tgrid(i) - ogrid);
   tol  = 1e-03;

   if any(tfac < tol)
     disp([' t = ' num2str(tgrid(i)) ' ' t_units])
     istor             = istor + 1;
     tarr(      istor) = tgrid(i);
     Taarr(     istor) = Ta;
     qaarr(     istor) = qa;
     Toarr(     istor) = To;
     qoarr(     istor) = qo;
     Tarr(    :,istor) = T';
     phi_iarr(:,istor) = phi_i';
     phi_uarr(:,istor) = phi_u';
     Karr(    :,istor) = K';
     Carr(    :,istor) = C';
     Jarr(    :,istor) = Jtilde';
     Jnetarr( :,istor) = Jnet';
     QSarr(   :,istor) = QS';
   end

% store values for next time step

   Tn  = T;
   qan = qa;
   qon = qo;
 end

% > Save variables to file

 Ofile = ['CVPMout/' experim '_cvpm.mat'];

% for a potential CVPM restart

 varout = 'planet site CS CS_limits rf R Dr dr varepR Lambda N Ar VP';

 varout = [varout ' t_units t_limits tstep ostep IEfac'];

 varout = [varout ' BCtypes BCfiles S_opt S0 hs'];

 varout = [varout ' Mtyp Km0 rhom cpm0 phi rho K C rhocp Psi1 Psi2 r1 r2 lambda solute xs0'];

 varout = [varout ' Mw phi_i phi_u phi_u_tab T_tab'];

 varout = [varout ' t Ta qa To qo T'];

% for CVPMview

 varout = [varout ' tarr Taarr qaarr Toarr qoarr Tarr'];

 varout = [varout ' phi_iarr phi_uarr Karr Carr'];
 
 varout = [varout ' Jarr Jnetarr QSarr'];

 eval(['save ' Ofile ' ' varout]);
