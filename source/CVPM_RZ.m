 function [] = CVPM_RZ(experim)

% Computational engine for the CVPM 2-D cylindrical case (RZ).
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
 innerBC_type = BCtypes{3};
 outerBC_type = BCtypes{4};

 upperBC_file = ['BCs/' BCfiles{1}];
 lowerBC_file = ['BCs/' BCfiles{2}];
 innerBC_file = ['BCs/' BCfiles{3}];
 outerBC_file = ['BCs/' BCfiles{4}];

 switch upperBC_type
 case 'T'
   [des,t_TsA,tunits_Ts,TsA,Tsmethod,R_TsA] = inputBC(upperBC_file,CS);
   t_TsA = convert_tunits(upperBC_file,t_TsA,tunits_Ts,t_units);
   TsR   = regridBC_2D(TsA,R_TsA,R,Tsmethod);	% Ts(t_TsA,R)
   Ts    = interp1(t_TsA,TsR,tmin,Tsmethod);    % Ts( tmin,R)
   qs    = zeros(size(Ts));
 case 'q'
   [des,t_qsA,tunits_qs,qsA,qsmethod,R_qsA] = inputBC(upperBC_file,CS);
   t_qsA = convert_tunits(upperBC_file,t_qsA,tunits_qs,t_units);
   qsR   = regridBC_2D(qsA,R_qsA,R,qsmethod);   % qs(t_qsA,R)
   qs    = interp1(t_qsA,qsR,tmin,qsmethod);    % qs( tmin,R)
   Ts    = zeros(size(qs));
 end

 switch lowerBC_type
 case 'T'
   [des,t_TbA,tunits_Tb,TbA,Tbmethod,R_TbA] = inputBC(lowerBC_file,CS);
   t_TbA = convert_tunits(lowerBC_file,t_TbA,tunits_Tb,t_units);
   TbR   = regridBC_2D(TbA,R_TbA,R,Tbmethod);   % Tb(t_TbA,R)
   Tb    = interp1(t_TbA,TbR,tmin,Tbmethod);    % Tb( tmin,R)
   qb    = zeros(size(Tb));
 case 'q'
   [des,t_qbA,tunits_qb,geoFA,qbmethod,R_qbA] = inputBC(lowerBC_file,CS);
   t_qbA = convert_tunits(lowerBC_file,t_qbA,tunits_qb,t_units);
   qbA   = -geoFA;
   qbR   = regridBC_2D(qbA,R_qbA,R,qbmethod);   % qb(t_qbA,R)
   qb    = interp1(t_qbA,qbR,tmin,qbmethod);    % qb( tmin,R)
   Tb    = zeros(size(qb));
 end

 switch innerBC_type
 case 'T'
   [des,t_TaA,tunits_Ta,dTaA,Tamethod,Z_TaA] = inputBC(innerBC_file,CS);
   t_TaA = convert_tunits(innerBC_file,t_TaA,tunits_Ta,t_units);
   dTaZ  = regridBC_2D(dTaA,Z_TaA,Z,Tamethod);  % dTa(t_TaA,Z)
   dTa   = interp1(t_TaA,dTaZ,tmin,Tamethod);   % dTa( tmin,Z)
   dTa   = dTa';
   qa    = zeros(size(dTa));
 case 'q'
   [des,t_qaA,tunits_qa,qaA,qamethod,Z_qaA] = inputBC(innerBC_file,CS);
   t_qaA = convert_tunits(innerBC_file,t_qaA,tunits_qa,t_units);
   qaZ   = regridBC_2D(qaA,Z_qaA,Z,qamethod);   % qa(t_qaA,Z)
   qa    = interp1(t_qaA,qaZ,tmin,qamethod);    % qa( tmin,Z)
   qa    = qa';
   dTa   = zeros(size(qa));
 case 'none'
   dTa   = zeros(size(Z));
   qa    = zeros(size(Z));
 end

 switch outerBC_type
 case 'T'
   [des,t_ToA,tunits_To,dToA,Tomethod,Z_ToA] = inputBC(outerBC_file,CS);
   t_ToA = convert_tunits(outerBC_file,t_ToA,tunits_To,t_units);
   dToZ  = regridBC_2D(dToA,Z_ToA,Z,Tomethod);  % dTo(t_ToA,Z)
   dTo   = interp1(t_ToA,dToZ,tmin,Tomethod);   % dTo( tmin,Z)
   dTo   = dTo';
   qo    = zeros(size(dTo));
 case 'q'
   [des,t_qoA,tunits_qo,qoA,qomethod,Z_qoA] = inputBC(outerBC_file,CS);
   t_qoA = convert_tunits(outerBC_file,t_qoA,tunits_qo,t_units);
   qoZ   = regridBC_2D(qoA,Z_qoA,Z,qomethod);   % qo(t_qoA,Z)
   qo    = interp1(t_qoA,qoZ,tmin,qomethod);    % qo( tmin,Z)
   qo    = qo';
   dTo   = zeros(size(qo));
 end

% setup computational time grid

 secs_tunit = dtunit2secs(t_units);     % seconds per t_unit

 dt    = secs_tunit * tstep;            % calculation interval (secs)
 Nt    = round((tmax-tmin)/tstep);      % number of time steps
 tgrid = tmin + tstep * (1:Nt);         % time grid (in t_units)

% setup output time grid

 if ostep < 0               % variable time step
   Ostep = abs(ostep);
   ogrid = [tmin    (Ostep:Ostep:10*Ostep)];
   Ostep = 10*Ostep;
   ogrid = [ogrid (2*Ostep:Ostep:10*Ostep)];
   Ostep = 10*Ostep;
   ogrid = [ogrid (2*Ostep:Ostep:10*Ostep)];
   Ostep = 10*Ostep;
   ogrid = [ogrid (2*Ostep:Ostep:1000*Ostep)];
   L     = ogrid <= tmax;
   ogrid = ogrid(L);
 else                       % uniform time step
   ogrid = (tmin:ostep:tmax);
 end
 Nout = length(ogrid);      % number of output times

% pre-allocate storage arrays

 tarr     = NaN*ones(      1,Nout);
 Tsarr    = NaN*ones(    N+1,Nout);
 qsarr    = NaN*ones(    N+1,Nout);
 Tbarr    = NaN*ones(    N+1,Nout);
 qbarr    = NaN*ones(    N+1,Nout);

 dTaarr   = NaN*ones(    M+1,Nout);
 qaarr    = NaN*ones(    M+1,Nout);
 dToarr   = NaN*ones(    M+1,Nout);
 qoarr    = NaN*ones(    M+1,Nout);

 Tarr     = NaN*ones(M+1,N+1,Nout);
 phi_iarr = NaN*ones(M+1,N+1,Nout);
 phi_uarr = NaN*ones(M+1,N+1,Nout);
 Karr     = NaN*ones(M+1,N+1,Nout);
 Carr     = NaN*ones(M+1,N+1,Nout);

% initial values for storage arrays

 istor                = 1;
 tarr(         istor) = tmin;
 Tsarr(      :,istor) = Ts;
 qsarr(      :,istor) = qs;
 Tbarr(      :,istor) = Tb;
 qbarr(      :,istor) = qb;

 dTaarr(     :,istor) = dTa;
 qaarr(      :,istor) = qa;
 dToarr(     :,istor) = dTo;
 qoarr(      :,istor) = qo;

 Tarr(     :,:,istor) = T;
 phi_iarr( :,:,istor) = phi_i;
 phi_uarr( :,:,istor) = phi_u;
 Karr(     :,:,istor) = K;
 Carr(     :,:,istor) = C;

% initialize additional arrays

 dphiudT = NaN*ones(size(phi_u));

% initial temperature values along the inner and outer boundaries

 Ta_init = T(:,  1,istor);
 To_init = T(:,N+1,istor);

% find the source-function integral (currently assumed to be independent of time)

 for j=1:N+1
   QS(:,j) = QSsub_Z(S_opt,S0(:,j),hs(:,j),zf);
 end

% ------ Begin Time Loop -------->

 Tn  = T;
 qsn = qs;
 qbn = qb;
 qan = qa;
 qon = qo;

 for i=1:Nt

   t = secs_tunit * tgrid(i);           % time (secs)

% find BCs

   switch upperBC_type
   case 'T'
     Ts = interp1(t_TsA,TsR,tgrid(i),Tsmethod);
     qs = zeros(size(Ts));
   case 'q'
     qs = interp1(t_qsA,qsR,tgrid(i),qsmethod);
     Ts = zeros(size(qs));
   end

   switch lowerBC_type
   case 'T'
     Tb = interp1(t_TbA,TbR,tgrid(i),Tbmethod);
     qb = zeros(size(Tb));
   case 'q'
     qb = interp1(t_qbA,qbR,tgrid(i),qbmethod);
     Tb = zeros(size(qb));
   end

   switch innerBC_type
   case 'T'
     dTa = interp1(t_TaA,dTaZ,tgrid(i),Tamethod);
     dTa = dTa';
     Ta  = Ta_init + dTa;
     qa  = zeros(size(dTa));
   case 'q'
     qa  = interp1(t_qaA,qaZ,tgrid(i),qamethod);
     qa  = qa';
     dTa = zeros(size(qa));
     Ta  = zeros(size(qa));
   case 'none'
     Ta  = Tn(:,1);
     qa  = zeros(size(Z));
   end

   switch outerBC_type
   case 'T'
     dTo = interp1(t_ToA,dToZ,tgrid(i),Tomethod);
     dTo = dTo';
     To  = To_init + dTo;
     qo  = zeros(size(dTo));
   case 'q'
     qo  = interp1(t_qoA,qoZ,tgrid(i),qomethod);
     qo  = qo';
     dTo = zeros(size(qo));
     To  = zeros(size(qo));
   end

   Tbc = {Ts, Tb, Ta, To};
   q   = {qs, qb, qa, qo};
   qn  = {qsn,qbn,qan,qon};

% find the volume fractions of ice and unfrozen water

   for k=1:M+1
     phi_u_tabS = squeeze(phi_u_tab(k,:,:));
     T_tabS     = squeeze(    T_tab(k,:,:));
     [phi_u(k,:),dphiudT(k,:)] ...
                = phiu_tableI(T(k,:),Mtyp(k,:),Mw(k,:),phi_u_tabS,T_tabS);
     phi_i(k,:) = phiiSub(Mtyp(k,:),Mw(k,:),phi_u(k,:));
   end

% find the thermal conductivity at the CV grid points

   for k=1:M+1
     K(k,:) = Ksub(T(k,:),Mtyp(k,:),Km0(k,:),phi(k,:),phi_i(k,:),phi_u(k,:),planet);
   end

% find effective conductivity at R- and Z-interfaces

   [KeR,KeZ] = Keff_RZ(K,varepR,varepZ);

% find the volumetric heat capacity

   for k=1:M+1
     [rhocp(k,:),C(k,:)] = Csub(T(k,:),Mtyp(k,:),rhom(k,:), ...
         cpm0(k,:),phi(k,:),phi_i(k,:),phi_u(k,:),dphiudT(k,:));
   end

% find discretization coefficients

   [aU,aUp,aD,aDp,aW,aWp,aE,aEp,aP,aPp,b] ...
      = DEcoefs_RZ(dt,f,rf,Lambda,Dz,Ar,Az,VP,C,KeR,KeZ,QS,BCtypes,q,qn);

% find new temperature field

   if mod(i,2)
     T = TDMA_RZh(aU,aP,aD,aW,aE,aUp,aPp,aDp,aWp,aEp,b,Tn,BCtypes,Tbc);
   else
     switch innerBC_type
     case {'T','q'}
       T = TDMA_RZh(aU,aP,aD,aW,aE,aUp,aPp,aDp,aWp,aEp,b,Tn,BCtypes,Tbc);
     case 'none'
       T = TDMA_RZv(aU,aP,aD,aW,aE,aUp,aPp,aDp,aWp,aEp,b,Tn,BCtypes,Tbc);
     end
   end

% estimate temperature on lower boundary if lowerBC_type = 'q'

   switch lowerBC_type
   case 'q'
     T(M+1,:) = T(M,:) - qb * dz(M+1) ./ KeZ(M+1,:);
     switch innerBC_type
     case 'T'
       T(M+1,1) = Ta(M+1);          % set corner pt back to Ta
     end
   end

% estimate temperature on inner boundary if innerBC_type = 'q'

   switch innerBC_type
   case 'q'
     T(:,1) = T(:,2) + qa * dr(2) ./ KeR(:,2);
   end

% estimate temperature on outer boundary if outerBC_type = 'q'

   switch outerBC_type
   case 'q'
     T(:,N+1) = T(:,N) - qo * dr(N+1) ./ KeR(:,N+1);
   end

% check numerical stability

   VPC           = VP .* C;
   AKeR          = Ar .* KeR;
   AKeZ          = Az .* KeZ;
   fac1          = NaN*ones(M+1,N+1);
   fac2          = NaN*ones(M+1,N+1);
   fac1(1:M,1:N) = AKeZ(1:M,1:N) + AKeZ(2:M+1,1:N);
   fac2(1:M,1:N) = AKeR(1:M,1:N) + AKeR(1:M,2:N+1);

   sfac = VPC ./ ((1-f) * (fac1 + fac2));

   switch innerBC_type
   case 'none'
     jmin = 2;
   case {'T','q'}
     jmin = 3;
   end
   minSfac = min(min(sfac(jmin:M,jmin:N)));

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
     istor                = istor + 1;
     tarr(         istor) = tgrid(i);
     Tsarr(      :,istor) = Ts;
     qsarr(      :,istor) = qs;
     Tbarr(      :,istor) = Tb;
     qbarr(      :,istor) = qb;
     dTaarr(     :,istor) = dTa;
     qaarr(      :,istor) = qa;
     dToarr(     :,istor) = dTo;
     qoarr(      :,istor) = qo;
     Tarr(     :,:,istor) = T;
     phi_iarr( :,:,istor) = phi_i;
     phi_uarr( :,:,istor) = phi_u;
     Karr(     :,:,istor) = K;
     Carr(     :,:,istor) = C;
   end

% store values for next time step

   Tn  = T;
   qsn = qs;
   qbn = qb;
   qan = qa;
   qon = qo;
 end

% > Save variables to file

 Ofile = ['CVPMout/' experim '_cvpm.mat'];

% for a potential CVPM restart

 varout = 'planet site CS CS_limits zf Z Dz dz varepZ M rf R Dr dr varepR Lambda N Ar Az VP';

 varout = [varout ' t_units t_limits tstep ostep IEfac'];

 varout = [varout ' BCtypes BCfiles S_opt S0 hs'];

 varout = [varout ' Mtyp Km0 rhom cpm0 phi rho K C rhocp solute xs0 lambda Psi r'];

 varout = [varout ' Mw phi_i phi_u phi_u_tab T_tab'];

 varout = [varout ' t Ts qs Tb qb Ta qa To qo T'];

% for CVPMview

 varout = [varout ' tarr Tsarr qsarr Tbarr qbarr dTaarr qaarr dToarr qoarr Tarr'];

 varout = [varout ' phi_iarr phi_uarr Karr Carr'];
 
 eval(['save ' Ofile ' ' varout]);
