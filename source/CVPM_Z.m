 function [] = CVPM_Z(experim)

% Computational engine for the CVPM 1-D depth case (Z).
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

 upperBC_type = BCtypes{1};
 lowerBC_type = BCtypes{2};

 upperBC_file = ['BCs/' BCfiles{1}];
 lowerBC_file = ['BCs/' BCfiles{2}];

 switch upperBC_type
 case 'T'
   [des,t_TsA,tunits_Ts,TsA,Tsmethod] = inputBC(upperBC_file,CS);
   t_TsA = convert_tunits(upperBC_file,t_TsA,tunits_Ts,t_units);
   Ts    = interp1(t_TsA,TsA,tmin,Tsmethod);
   qs    = 0;
 case 'q'
   [des,t_qsA,tunits_qs,qsA,qsmethod] = inputBC(upperBC_file,CS);
   t_qsA = convert_tunits(upperBC_file,t_qsA,tunits_qs,t_units);
   qs    = interp1(t_qsA,qsA,tmin,qsmethod);
   Ts    = 0;
 end

 switch lowerBC_type
 case 'T'
   [des,t_TbA,tunits_Tb,TbA,Tbmethod] = inputBC(lowerBC_file,CS);
   t_TbA = convert_tunits(lowerBC_file,t_TbA,tunits_Tb,t_units);
   Tb    = interp1(t_TbA,TbA,tmin,Tbmethod);
   qb    = 0;
 case 'q'
   [des,t_qbA,tunits_qb,geoFA,qbmethod] = inputBC(lowerBC_file,CS);
   t_qbA = convert_tunits(lowerBC_file,t_qbA,tunits_qb,t_units);
   qbA   = -geoFA;
   qb    = interp1(t_qbA,qbA,tmin,qbmethod);
   Tb    = 0;
 end

% setup time grid

 secs_tunit = dtunit2secs(t_units);     % seconds per t_unit

 dt      = secs_tunit * tstep;          % calculation interval (secs)
 Nt      = round((tmax-tmin)/tstep);    % number of time steps
 tgrid   = tmin + tstep * (1:Nt);       % time grid (in t_units)
 outFreq = round(ostep/tstep);          % output frequency
 Nout    = (Nt / outFreq) + 1;          % number of output times

% pre-allocate storage arrays

 tarr     = NaN*ones(  1,Nout);
 Tsarr    = NaN*ones(  1,Nout);
 qsarr    = NaN*ones(  1,Nout);
 Tbarr    = NaN*ones(  1,Nout);
 qbarr    = NaN*ones(  1,Nout);
 Tarr     = NaN*ones(M+1,Nout);
 phi_iarr = NaN*ones(M+1,Nout);
 phi_uarr = NaN*ones(M+1,Nout);
 Karr     = NaN*ones(M+1,Nout);
 Carr     = NaN*ones(M+1,Nout);
 Jarr     = NaN*ones(M+1,Nout);
 Jnetarr  = NaN*ones(M+1,Nout);
 QSarr    = NaN*ones(M+1,Nout);

% initial values for storage arrays

 Ke       = Keff(K,varepZ);
 [J,Jnet] = Jsub_Z(BCtypes,T,Ke,dz,qb);

 istor             = 1;
 tarr(      istor) = tmin;
 Tsarr(     istor) = Ts;
 qsarr(     istor) = qs;
 Tbarr(     istor) = Tb;
 qbarr(     istor) = qb;
 Tarr(    :,istor) = T;
 phi_iarr(:,istor) = phi_i;
 phi_uarr(:,istor) = phi_u;
 Karr(    :,istor) = K;
 Carr(    :,istor) = C;
 Jarr(    :,istor) = J;
 Jnetarr( :,istor) = Jnet;

% find the source-function integral (currently assumed to be independent of time)

 QS = QSsub_Z(S_opt,S0,hs,zf);

% ------ Begin Time Loop -------->

 Tn  = T;
 qsn = qs;
 qbn = qb;

 for i=1:Nt

   t = secs_tunit * tgrid(i);           % time (secs)

% find BCs

   switch upperBC_type
   case 'T'
     Ts = interp1(t_TsA,TsA,tgrid(i),Tsmethod);
     qs = 0;
   otherwise
     Ts = 0;
     qs = interp1(t_qsA,qsA,tgrid(i),qsmethod);
   end

   switch lowerBC_type
   case 'T'
     Tb = interp1(t_TbA,TbA,tgrid(i),Tbmethod);
     qb = 0;
   case 'q';
     Tb = 0;
     qb = interp1(t_qbA,qbA,tgrid(i),qbmethod);
   end

   Tbc = [Ts  Tb];
   q   = [qs  qb];
   qn  = [qsn qbn];

% find the volume fractions of ice and unfrozen water

   [phi_u,dphiudT] = phiu_tableI(T,Mtyp,Mw,phi_u_tab,T_tab);
   phi_i           = phiiSub(Mtyp,Mw,phi_u);

% find the thermal conductivity at the CV grid points

   K = Ksub(T,Mtyp,Km0,phi,phi_i,phi_u,planet);

% find effective thermal conductivity at CV interfaces

   Ke = Keff(K,varepZ);

% find the volumetric heat capacity

   [rhocp,C] = Csub(T,Mtyp,rhom,cpm0,phi,phi_i,phi_u,dphiudT);

% find discretization coefficients

   [aU,aUp,aD,aDp,aP,aPp,b] = DEcoefs_Z(dt,f,Az,VP,C,Ke,QS,BCtypes,q,qn);

% find new temperature field

   T = TDMA_Z(aU,aP,aD,aUp,aPp,aDp,b,Tn,BCtypes,Tbc);

% estimate temperature on lower boundary if lowerBC_type = 'q'

   switch lowerBC_type
   case 'q'
     T(M+1) = T(M) - qb * dz(M+1) / Ke(M+1);
   end

% diagnostics: find diffusive heat-flux across the interfaces

   [J,Jnet] = Jsub_Z(BCtypes,T,Ke,dz,qb);

% check numerical stability

   VPC      = VP .* C;
   AKeZ     = Az .* Ke;
   fac      = NaN*ones(M+1,1);
   fac(1:M) = AKeZ(1:M) + AKeZ(2:M+1);

   sfac     = VPC ./ ((1-f) * fac);
   minSfac  = min(sfac(2:M));

   if minSfac < dt
     disp(' ')
     disp(['Warning: time step is too big for numerical stability, t = ' num2str(tgrid(i)) ' (' t_units ')'])
     secs_tunit = dtunit2secs(t_units);     % seconds per t_unit
     minSfac    = minSfac / secs_tunit;     % convert to t_units
     disp(' ')
     disp(['dt should be less than ' num2str(minSfac) ' (' t_units ') to guarantee numerical stability.'])
   end

% store results for output

   tratio = i/outFreq;

   if isequal(tratio,floor(tratio))
     if tratio < 4 || i == Nt
       disp([' t = ' num2str(tgrid(i)) ' ' t_units])
     end
     istor             = istor + 1;
     tarr(      istor) = tgrid(i);
     Tsarr(     istor) = Ts;
     qsarr(     istor) = qs;
     Tbarr(     istor) = Tb;
     qbarr(     istor) = qb;
     Tarr(    :,istor) = T;
     phi_iarr(:,istor) = phi_i;
     phi_uarr(:,istor) = phi_u;
     Karr(    :,istor) = K;
     Carr(    :,istor) = C;
     Jarr(    :,istor) = J;
     Jnetarr( :,istor) = Jnet;
     QSarr(   :,istor) = QS;
   end

% store values for next time step

   Tn  = T;
   qsn = qs;
   qbn = qb;
 end

% > Save variables to file

 Ofile = ['CVPMout/' experim '_cvpm.mat'];

% for a potential CVPM restart

 varout = 'planet site CS CS_limits zf Z Dz dz varepZ M Az VP';

 varout = [varout ' t_units t_limits tstep ostep IEfac'];

 varout = [varout ' BCtypes BCfiles S_opt S0 hs'];

 varout = [varout ' Mtyp Km0 rhom cpm0 phi rho K C rhocp Psi1 Psi2 r1 r2 lambda solute xs0'];

 varout = [varout ' Mw phi_i phi_u phi_u_tab T_tab'];

 varout = [varout ' t Ts qs Tb qb T'];

% for CVPMview

 varout = [varout ' tarr Tsarr qsarr Tbarr qbarr Tarr'];

 varout = [varout ' phi_iarr phi_uarr Karr Carr'];
 
 varout = [varout ' Jarr Jnetarr QSarr'];

 eval(['save ' Ofile ' ' varout]);
