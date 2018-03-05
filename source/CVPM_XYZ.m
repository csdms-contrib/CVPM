 function [] = CVPM_XYZ(experim)

% Computational engine for the CVPM 3-D cartesian case (XYZ).
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
 xleftBC_type = BCtypes{3};
 xrighBC_type = BCtypes{4};
 yleftBC_type = BCtypes{5};
 yrighBC_type = BCtypes{6};

 upperBC_file = ['BCs/' BCfiles{1}];
 lowerBC_file = ['BCs/' BCfiles{2}];
 xleftBC_file = ['BCs/' BCfiles{3}];
 xrighBC_file = ['BCs/' BCfiles{4}];
 yleftBC_file = ['BCs/' BCfiles{5}];
 yrighBC_file = ['BCs/' BCfiles{6}];

 switch upperBC_type
 case 'T'
   [des,t_TsA,tunits_Ts,TsA,Tsmethod,X_TsA,Y_TsA] = inputBC(upperBC_file,CS);
   t_TsA = convert_tunits(upperBC_file,t_TsA,tunits_Ts,t_units);
   TsXY  = regridBC_3D(TsA,X_TsA,Y_TsA,X,Y,Tsmethod);	% Ts(t_TsA,X,Y)
   Ts    = interpBC_3D(t_TsA,TsXY,tmin,Tsmethod);       % Ts( tmin,X,Y)
   qs    = zeros(size(Ts));
 case 'q'
   [des,t_qsA,tunits_qs,qsA,qsmethod,X_qsA,Y_qsA] = inputBC(upperBC_file,CS);
   t_qsA = convert_tunits(upperBC_file,t_qsA,tunits_qs,t_units);
   qsXY  = regridBC_3D(qsA,X_qsA,Y_qsA,X,Y,qsmethod);   % qs(t_qsA,X,Y)
   qs    = interpBC_3D(t_qsA,qsXY,tmin,qsmethod);       % qs( tmin,X,Y)
   Ts    = zeros(size(qs));
 end

 switch lowerBC_type
 case 'T'
   [des,t_TbA,tunits_Tb,TbA,Tbmethod,X_TbA,Y_TbA] = inputBC(lowerBC_file,CS);
   t_TbA = convert_tunits(lowerBC_file,t_TbA,tunits_Tb,t_units);
   TbXY  = regridBC_3D(TbA,X_TbA,Y_TbA,X,Y,Tbmethod);	% Tb(t_TbA,X,Y)
   Tb    = interpBC_3D(t_TbA,TbXY,tmin,Tbmethod);       % Tb( tmin,X,Y)
   qb    = zeros(size(Tb));
 case 'q'
   [des,t_qbA,tunits_qb,geoFA,qbmethod,X_qbA,Y_qbA] = inputBC(lowerBC_file,CS);
   t_qbA = convert_tunits(lowerBC_file,t_qbA,tunits_qb,t_units);
   qbA   = -geoFA;
   qbXY  = regridBC_3D(qbA,X_qbA,Y_qbA,X,Y,qbmethod);	% qb(t_qbA,X,Y)
   qb    = interpBC_3D(t_qbA,qbXY,tmin,qbmethod);       % qb( tmin,X,Y)
   Tb    = zeros(size(qb));
 end

 switch xleftBC_type
 case 'T'
   [des,t_TaA,tunits_Ta,dTaA,Tamethod,Y_TaA,Z_TaA] = inputBC(xleftBC_file,CS);
   t_TaA = convert_tunits(xleftBC_file,t_TaA,tunits_Ta,t_units);
   dTaYZ = regridBC_3D(dTaA,Y_TaA,Z_TaA,Y,Z,Tamethod);	% dTa(t_TaA,Y,Z)
   dTa   = interpBC_3D(t_TaA,dTaYZ,tmin,Tamethod);      % dTa( tmin,Y,Z)
   dTa   = dTa';
   qa    = zeros(size(dTa));
 case 'q'
   [des,t_qaA,tunits_qa,qaA,qamethod,Y_qaA,Z_qaA] = inputBC(xleftBC_file,CS);
   t_qaA = convert_tunits(xleftBC_file,t_qaA,tunits_qa,t_units);
   qaYZ  = regridBC_3D(qaA,Y_qaA,Z_qaA,Y,Z,qamethod);	% qa(t_qaA,Y,Z)
   qa    = interpBC_3D(t_qaA,qaYZ,tmin,qamethod);       % qa( tmin,Z)
   qa    = qa';
   dTa   = zeros(size(qa));
 end

 switch xrighBC_type
 case 'T'
   [des,t_ToA,tunits_To,dToA,Tomethod,Y_ToA,Z_ToA] = inputBC(xrighBC_file,CS);
   t_ToA = convert_tunits(xrighBC_file,t_ToA,tunits_To,t_units);
   dToYZ = regridBC_3D(dToA,Y_ToA,Z_ToA,Y,Z,Tomethod);	% dTo(t_ToA,Y,Z)
   dTo   = interpBC_3D(t_ToA,dToYZ,tmin,Tomethod);      % dTo( tmin,Z)
   dTo   = dTo';
   qo    = zeros(size(dTo));
 case 'q'
   [des,t_qoA,tunits_qo,qoA,qomethod,Y_qoA,Z_qoA] = inputBC(xrighBC_file,CS);
   t_qoA = convert_tunits(xrighBC_file,t_qoA,tunits_qo,t_units);
   qoYZ  = regridBC_3D(qoA,Y_qoA,Z_qoA,Y,Z,qomethod);	% qo(t_qoA,Y,Z)
   qo    = interpBC_3D(t_qoA,qoYZ,tmin,qomethod);       % qo( tmin,Z)
   qo    = qo';
   dTo   = zeros(size(qo));
 end

 switch yleftBC_type
 case 'T'
   [des,t_TcA,tunits_Tc,dTcA,Tcmethod,X_TcA,Z_TcA] = inputBC(yleftBC_file,CS);
   t_TcA = convert_tunits(yleftBC_file,t_TcA,tunits_Tc,t_units);
   dTcXZ = regridBC_3D(dTcA,X_TcA,Z_TcA,X,Z,Tcmethod);	% dTc(t_TcA,X,Z)
   dTc   = interpBC_3D(t_TcA,dTcXZ,tmin,Tcmethod);      % dTc( tmin,X,Z)
   dTc   = dTc';
   qc    = zeros(size(dTc));
 case 'q'
   [des,t_qcA,tunits_qc,qcA,qcmethod,X_qcA,Z_qcA] = inputBC(yleftBC_file,CS);
   t_qcA = convert_tunits(yleftBC_file,t_qcA,tunits_qc,t_units);
   qcXZ = regridBC_3D(qcA,X_qcA,Z_qcA,X,Z,qcmethod);	% qc(t_qcA,X,Z)
   qc   = interpBC_3D(t_qcA,qcXZ,tmin,qcmethod);        % qc( tmin,Z)
   qc   = qc';
   dTc  = zeros(size(qc));
 end

 switch yrighBC_type
 case 'T'
   [des,t_TdA,tunits_Td,dTdA,Tdmethod,X_TdA,Z_TdA] = inputBC(yrighBC_file,CS);
   t_TdA = convert_tunits(yrighBC_file,t_TdA,tunits_Td,t_units);
   dTdXZ = regridBC_3D(dTdA,X_TdA,Z_TdA,X,Z,Tdmethod);	% dTd(t_TdA,X,Z)
   dTd   = interpBC_3D(t_TdA,dTdXZ,tmin,Tdmethod);      % dTd( tmin,Z)
   dTd   = dTd';
   qd    = zeros(size(dTd));
 case 'q'
   [des,t_qdA,tunits_qd,qdA,qdmethod,X_qdA,Z_qdA] = inputBC(yrighBC_file,CS);
   t_qdA = convert_tunits(yrighBC_file,t_qdA,tunits_qd,t_units);
   qdXZ = regridBC_3D(qdA,X_qdA,Z_qdA,X,Z,qdmethod);	% qd(t_qdA,X,Z)
   qd   = interpBC_3D(t_qdA,qdXZ,tmin,qdmethod);        % qd( tmin,Z)
   qd   = qd';
   dTd  = zeros(size(qd));
 end

% setup time grid

 secs_tunit = dtunit2secs(t_units);     % seconds per t_unit

 dt      = secs_tunit * tstep;          % calculation interval (secs)
 Nt      = round((tmax-tmin)/tstep);    % number of time steps
 tgrid   = tmin + tstep * (1:Nt);       % time grid (in t_units)

 ogrid = (tmin:ostep:tmax);             % uniform time step
 Nout  = length(ogrid);                 % number of output times

% pre-allocate storage arrays

 tarr     = NaN*ones(      1,Nout);
 Tsarr    = NaN*ones(N+1,L+1,Nout);
 qsarr    = NaN*ones(N+1,L+1,Nout);
 Tbarr    = NaN*ones(N+1,L+1,Nout);
 qbarr    = NaN*ones(N+1,L+1,Nout);

 dTaarr   = NaN*ones(M+1,L+1,Nout);
 qaarr    = NaN*ones(M+1,L+1,Nout);
 dToarr   = NaN*ones(M+1,L+1,Nout);
 qoarr    = NaN*ones(M+1,L+1,Nout);

 Tarr     = NaN*ones(M+1,N+1,L+1,Nout);
 phi_iarr = NaN*ones(M+1,N+1,L+1,Nout);
 phi_uarr = NaN*ones(M+1,N+1,L+1,Nout);
 Karr     = NaN*ones(M+1,N+1,L+1,Nout);
 Carr     = NaN*ones(M+1,N+1,L+1,Nout);

% initial values for storage arrays

 istor                  = 1;
 tarr(           istor) = tmin;
 Tsarr(      :,:,istor) = Ts;
 qsarr(      :,:,istor) = qs;
 Tbarr(      :,:,istor) = Tb;
 qbarr(      :,:,istor) = qb;

 dTaarr(     :,:,istor) = dTa;
 qaarr(      :,:,istor) = qa;
 dToarr(     :,:,istor) = dTo;
 qoarr(      :,:,istor) = qo;

 Tarr(     :,:,:,istor) = T;
 phi_iarr( :,:,:,istor) = phi_i;
 phi_uarr( :,:,:,istor) = phi_u;
 Karr(     :,:,:,istor) = K;
 Carr(     :,:,:,istor) = C;

% initialize additional arrays

 dphiudT = NaN*ones(size(phi_u));

% initial temperature values along the X,Y boundaries

 Ta_init = NaN*ones(size(dTa));
 To_init = NaN*ones(size(dTo));
 Tc_init = NaN*ones(size(dTc));
 Td_init = NaN*ones(size(dTd));

 Ta_init(:,:) = T(:,  1,:,istor);
 To_init(:,:) = T(:,N+1,:,istor);
 Tc_init(:,:) = T(:,:,  1,istor);
 Td_init(:,:) = T(:,:,L+1,istor);

% find the source-function integral (currently assumed to be independent of time)

 for j=1:N+1
   for jy=1:L+1
     QS(:,j,jy) = QSsub_Z(S_opt,S0(:,j,jy),hs(:,j,jy),zf);
   end
 end

% ------ Begin Time Loop -------->

 Tn  = T;
 qsn = qs;
 qbn = qb;
 qan = qa;
 qon = qo;
 qcn = qc;
 qdn = qd;

 for i=1:Nt

   t = secs_tunit * tgrid(i);           % time (secs)

% find BCs

   switch upperBC_type
   case 'T'
     Ts = interpBC_3D(t_TsA,TsXY,tgrid(i),Tsmethod);
     qs = zeros(size(Ts));
   case 'q'
     qs = interpBC_3D(t_qsA,qsXY,tgrid(i),qsmethod);
     Ts = zeros(size(qs));
   end

   switch lowerBC_type
   case 'T'
     Tb = interpBC_3D(t_TbA,TbXY,tgrid(i),Tbmethod);
     qb = zeros(size(Tb));
   case 'q'
     qb = interpBC_3D(t_qbA,qbXY,tgrid(i),qbmethod);
     Tb = zeros(size(qb));
   end

   switch xleftBC_type
   case 'T'
     dTa = interpBC_3D(t_TaA,dTaYZ,tgrid(i),Tamethod);
     dTa = dTa';
     Ta  = Ta_init + dTa;
     qa  = zeros(size(dTa));
   case 'q'
     qa  = interpBC_3D(t_qaA,qaYZ,tgrid(i),qamethod);
     qa  = qa';
     dTa = zeros(size(qa));
     Ta  = zeros(size(qa));
   end

   switch xrighBC_type
   case 'T'
     dTo = interpBC_3D(t_ToA,dToYZ,tgrid(i),Tomethod);
     dTo = dTo';
     To  = To_init + dTo;
     qo  = zeros(size(dTo));
   case 'q'
     qo  = interpBC_3D(t_qoA,qoYZ,tgrid(i),qomethod);
     qo  = qo';
     dTo = zeros(size(qo));
     To  = zeros(size(qo));
   end

   switch yleftBC_type
   case 'T'
     dTc = interpBC_3D(t_TcA,dTcXZ,tgrid(i),Tcmethod);
     dTc = dTc';
     Tc  = Tc_init + dTc;
     qc  = zeros(size(dTc));
   case 'q'
     qc  = interpBC_3D(t_qcA,qcXZ,tgrid(i),qcmethod);
     qc  = qc';
     dTc = zeros(size(qc));
     Tc  = zeros(size(qc));
   end

   switch yrighBC_type
   case 'T'
     dTd = interpBC_3D(t_TdA,dTdXZ,tgrid(i),Tdmethod);
     dTd = dTd';
     Td  = Td_init + dTd;
     qd  = zeros(size(dTd));
   case 'q'
     qd  = interpBC_3D(t_qdA,qdXZ,tgrid(i),qdmethod);
     qd  = qd';
     dTd = zeros(size(qd));
     Td  = zeros(size(qd));
   end

   Tbc = {Ts, Tb, Ta, To, Tc, Td};
   q   = {qs, qb, qa, qo, qc, qd};
   qn  = {qsn,qbn,qan,qon,qcn,qdn};

% find the volume fractions of ice and unfrozen water

   for jy=1:L+1
     for k=1:M+1
       phi_u_tabS    = squeeze(phi_u_tab(k,:,jy,:));
       T_tabS        = squeeze(    T_tab(k,:,jy,:));
       [phi_u(k,:,jy),dphiudT(k,:,jy)] ...
                     = phiu_tableI(T(k,:,jy),Mtyp(k,:,jy),Mw(k,:,jy),phi_u_tabS,T_tabS);
       phi_i(k,:,jy) = phiiSub(Mtyp(k,:,jy),Mw(k,:,jy),phi_u(k,:,jy));
     end
   end

% find the thermal conductivity at the CV grid points

   for jy=1:L+1
     for k=1:M+1
       K(k,:,jy) = Ksub(T(k,:,jy),Mtyp(k,:,jy),Km0(k,:,jy),phi(k,:,jy), ...
                                 phi_i(k,:,jy),phi_u(k,:,jy),planet);
     end
   end

% effective conductivity at X-, Y-, Z-interfaces

   [KeX,KeY,KeZ] = Keff_XYZ(K,varepX,varepY,varepZ);

% find the volumetric heat capacity

   for jy=1:L+1
     for k=1:M+1
       [rhocp(k,:,jy),C(k,:,jy)] = Csub(T(k,:,jy),Mtyp(k,:,jy),rhom(k,:,jy), ...
           cpm0(k,:,jy),phi(k,:,jy),phi_i(k,:,jy),phi_u(k,:,jy),dphiudT(k,:,jy));
     end
   end

% find discretization coefficients

   [aU,aUp,aD,aDp,aW,aWp,aE,aEp,aS,aSp,aN,aNp,aP,aPp,b] ...
      = DEcoefs_XYZ(dt,f,Dx,Dy,Dz,Ax,Ay,Az,VP,C,KeX,KeY,KeZ,QS,BCtypes,q,qn);

% find new temperature field

   T = TDMA_XYZ(aU,aP,aD,aW,aE,aS,aN,aUp,aPp,aDp,aWp,aEp,aSp,aNp,b,Tn, ...
         BCtypes,Tbc);

% estimate temperature on lower boundary if lowerBC_type = 'q'

   switch lowerBC_type
   case 'q'
     T(M+1,:,:) = T(M,:,:) - permute(qb,[3 1 2]) * dz(M+1) ./ KeZ(M+1,:,:);
   end

% estimate temperature on inner boundary if xleftBC_type = 'q'

   switch xleftBC_type
   case 'q'
     T(:,1,:) = T(:,2,:) + permute(qa,[1 3 2]) * dx(2) ./ KeX(:,2,:);
   end

% estimate temperature on outer boundary if xrighBC_type = 'q'

   switch xrighBC_type
   case 'q'
     T(:,N+1,:) = T(:,N,:) - permute(qo,[1 3 2]) * dx(N+1) ./ KeX(:,N+1,:);
   end

% estimate temperature on inner boundary if yleftBC_type = 'q'

   switch yleftBC_type
   case 'q'
     T(:,:,1) = T(:,:,2) + permute(qc,[1 2 3]) * dy(2) ./ KeY(:,:,2);
   end

% estimate temperature on outer boundary if yrighBC_type = 'q'

   switch yrighBC_type
   case 'q'
     T(:,:,L+1) = T(:,:,L) - permute(qd,[1 2 3]) * dy(L+1) ./ KeY(:,:,L+1);
   end

% check numerical stability

   VPC  = VP .* C;
   AKeX = Ax .* KeX;
   AKeY = Ay .* KeY;
   AKeZ = Az .* KeZ;
   fac1 = NaN*ones(M+1,N+1,L+1);
   fac2 = NaN*ones(M+1,N+1,L+1);
   fac3 = NaN*ones(M+1,N+1,L+1);
   fac1(1:M,1:N,1:L) = AKeZ(1:M,1:N,1:L) + AKeZ(2:M+1,1:N,1:L);
   fac2(1:M,1:N,1:L) = AKeX(1:M,1:N,1:L) + AKeX(1:M,2:N+1,1:L);
   fac3(1:M,1:N,1:L) = AKeY(1:M,1:N,1:L) + AKeY(1:M,1:N,2:L+1);

   sfac    = VPC ./ ((1-f) * (fac1 + fac2 + fac3));
   minSfac = min(min(min(sfac(2:M,2:N,2:L))));

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
     istor                 = istor + 1;
     tarr(          istor) = tgrid(i);
     Tsarr(     :,:,istor) = Ts;
     qsarr(     :,:,istor) = qs;
     Tbarr(     :,:,istor) = Tb;
     qbarr(     :,:,istor) = qb;

     dTaarr(    :,:,istor) = dTa;
     qaarr(     :,:,istor) = qa;
     dToarr(    :,:,istor) = dTo;
     qoarr(     :,:,istor) = qo;

     Tarr(    :,:,:,istor) = T;
     phi_iarr(:,:,:,istor) = phi_i;
     phi_uarr(:,:,:,istor) = phi_u;
     Karr(    :,:,:,istor) = K;
     Carr(    :,:,:,istor) = C;
   end

% store values for next time step

   Tn  = T;
   qsn = qs;
   qbn = qb;
   qan = qa;
   qon = qo;
   qcn = qc;
   qdn = qd;
 end

% > Save variables to file

 Ofile = ['CVPMout/' experim '_cvpm.mat'];

% for a potential CVPM restart

 varout = 'planet site CS CS_limits zf Z Dz dz varepZ M xf X Dx dx varepX N yf Y Dy dy varepY L Ax Ay Az VP';

 varout = [varout ' t_units t_limits tstep ostep IEfac'];

 varout = [varout ' BCtypes BCfiles S_opt S0 hs'];

 varout = [varout ' Mtyp Km0 rhom cpm0 phi rho K C rhocp Psi1 Psi2 r1 r2 lambda solute xs0'];

 varout = [varout ' Mw phi_i phi_u phi_u_tab T_tab'];

 varout = [varout ' t Ts qs Tb qb Ta qa To qo Tc qc Td qd T'];

% for CVPMview

 varout = [varout ' tarr Tsarr qsarr Tbarr qbarr dTaarr qaarr dToarr qoarr Tarr'];

 varout = [varout ' phi_iarr phi_uarr Karr Carr'];
 
 eval(['save ' Ofile ' ' varout]);
