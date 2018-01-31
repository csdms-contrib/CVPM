 function [] = CPS_RZ(experim,opt)

% CVPM pre-processor for the 2-D cylindrical case (RZ).

% Indices:
%   i   time
%   j   R-index
%   k   Z-index

% Notes:
%   (1) For 'local' scale problems, the initial temperature field is 
%       calculated using the material properties found near the edge of
%       of the domain (assuming the initial field is not imported from 
%       from an IC_file).  For 'regional' scale problems, the initial
%       temperature field is found using material properties found at 
%       each location in the grid.
% ______________________________________________

% Written by:

%   Gary Clow
%   Institute of Arctic and Alpine Research
%   University of Colorado
%   Boulder, Colorado USA
%   Email: gary.clow@colorado.edu
% ______________________________________________

 disp(' ')
 disp('* Running CPS *')

 pos = set_screen(0);

% diagnostics

 dflag = 1;         % (0) quiet mode, (1) dumps various fields to screen

% parameters

 Tf0    = 0.01;     % triple point temperature (C)
 Ttol   = 1e-14;    % tolerance for initial temperature field
 rhotol = 1e-10;    % tolerance for density field

% input namelist defining the experiment

 Nfile = ['namelists/' experim '.namelist'];

 [planet,site,CS,Pscale,Bdepth,CS_limits,t_units,t_limits,tstep,ostep,initT_opt, ...
    ICfile,BCtypes,BCfiles,S_opt,C_opt,P_opt,solute,IEfac] = input_namelist_RZ(Nfile);

% set limits and BC/IC file names

 switch opt(1)
 case 1
   Zlayer_file = ['geo/' site '_Zlayers.txt'];
 case 2
   Zlayer_file = ['geo/' site '_Zlayers.mat'];
 case 3
   Zlayer_file = ['tmp/' site '_Zlayers_tmp.mat'];
 end
 Rlayer_file = ['geo/' site '_Rlayers.txt'];

 tmin = t_limits(1);
 Zmin = CS_limits(1);
 Zmax = CS_limits(2);
 Rmin = CS_limits(3);
 Rmax = CS_limits(4);

 upperBC_type = BCtypes{1};
 lowerBC_type = BCtypes{2};
 innerBC_type = BCtypes{3};
 outerBC_type = BCtypes{4};

 IC_file      = ['ICs/' ICfile];
 upperBC_file = ['BCs/' BCfiles{1}];
 lowerBC_file = ['BCs/' BCfiles{2}];
 innerBC_file = ['BCs/' BCfiles{3}];
 outerBC_file = ['BCs/' BCfiles{4}];

 disp(' ')
 disp(['Zlayer_file = ' Zlayer_file])
 disp(['Rlayer_file = ' Rlayer_file])
 disp(' ')
 disp(['upperBC_type = ' upperBC_type])
 disp(['lowerBC_type = ' lowerBC_type])
 disp(['innerBC_type = ' innerBC_type])
 disp(['outerBC_type = ' outerBC_type])

% > Setup spatial grid ---------------------------------

% input layers

 [zfuL,zfdL,dzL,MtypZL,MarrZL] = input_layers(Zlayer_file);
 [rfwL,rfeL,drL,MtypRL,MarrRL] = input_layers(Rlayer_file);

% setup CV grid

 [zf,Z,Dz,dz,varepZ,MtypZ,MarrZ,M] = Zgrid(zfuL,zfdL,dzL,MtypZL,MarrZL,Zmin,Zmax);
 [rf,R,Dr,dr,varepR,Lambda,MtypR,MarrR,N] = Rgrid(rfwL,rfeL,drL,MtypRL,MarrRL, ...
        innerBC_type,Rmin,Rmax);

% define geometric factors

 VP = repmat(Lambda,[M+1 1]) .* repmat(Dz,[1 N+1]);
 Az = repmat(Lambda,[M+1 1]) ./ repmat(dz,[1 N+1]);
 Ar = repmat(rf,[M+1 1]) .* repmat(Dz,[1 N+1]) ./ repmat(dr,[M+1 1]);

 switch innerBC_type
 case 'T'
   Ar(:,2) = Ar(:,2) + 0.5;
 end
 switch outerBC_type
 case 'T'
   Ar(:,N+1) = Ar(:,N+1) - 0.5;
 end

% initialize material properties

 Mtyp   = NaN*ones(M+1,N+1);
 Km0    = NaN*ones(M+1,N+1);
 rhom   = NaN*ones(M+1,N+1);
 cpm0   = NaN*ones(M+1,N+1);
 S0     = NaN*ones(M+1,N+1);
 hs     = NaN*ones(M+1,N+1);
 QS     = NaN*ones(M+1,N+1);
 phi    =    zeros(M+1,N+1);
 Sr     =     ones(M+1,N+1);
 xs0    =    zeros(M+1,N+1);
 lambda = NaN*ones(M+1,N+1);
 Psi    =     ones(M+1,N+1);
 r      = NaN*ones(M+1,N+1);

 for j=1:N+1
   switch MtypR(j)
   case 99              % use what's specified in the *_Zlayers file
     Mtyp(  :,j) = MtypZ;
     Km0(   :,j) = init_Km(         MtypZ,MarrZ, 1,'Z');    % matrix conductivity at 0 C
     rhom(  :,j) = init_rhom(       MtypZ,MarrZ, 2,'Z');    % matrix density
     cpm0(  :,j) = init_cpm(        MtypZ,MarrZ, 3,'Z');    % matrix specific heat at 20 C
     S0(    :,j) = init_S0(               MarrZ, 4,'Z');    % heat-source extrapolated to surface
     hs(    :,j) = init_hs(               MarrZ, 5,'Z');    % heat-source length scale
     phi(   :,j) = init_phi(C_opt,Z,MtypZ,MarrZ, 6,'Z');    % porosity
     Sr(    :,j) = init_Sr(         MtypZ,MarrZ, 9,'Z');    % degree of pore saturation
     xs0(   :,j) = init_xs0(        MtypZ,MarrZ,10,'Z');    % mole fraction of solutes (with no ice)
     lambda(:,j) = init_lambda(     MtypZ,MarrZ,11,'Z');    % interfacial melting paramater
     Psi(   :,j) = init_Psi(        MtypZ,MarrZ,12,'Z');    % volume fraction of pores associated with particles of radius r
     r(     :,j) = init_r(          MtypZ,MarrZ,13,'Z');    % effective radius of matrix particles

   otherwise        
     L = Z <= Bdepth;   % use what's in *_Rlayers file above borehole bottom
     Mtyp(  L,j) = MtypR(j);
     Km0(   L,j) = MarrR( 1,j);
     rhom(  L,j) = MarrR( 2,j);
     cpm0(  L,j) = MarrR( 3,j);
     S0(    L,j) = MarrR( 4,j);
     hs(    L,j) = MarrR( 5,j);
     phi(   L,j) = MarrR( 6,j);
     Sr(    L,j) = MarrR( 9,j);
     xs0(   L,j) = MarrR(10,j);
     lambda(L,j) = MarrR(11,j);
     Psi(   L,j) = MarrR(12,j);
     r(     L,j) = MarrR(13,j);

     L = Z > Bdepth;   % use what's in *_Zlayers file below borehole bottom
     if any(L)
       Mtyp(  L,j) = MtypZ(L); 
       Km0(   L,j) = init_Km(            MtypZ(L),MarrZ(L,:), 1,'Z');
       rhom(  L,j) = init_rhom(          MtypZ(L),MarrZ(L,:), 2,'Z');
       cpm0(  L,j) = init_cpm(           MtypZ(L),MarrZ(L,:), 3,'Z');
       S0(    L,j) = init_S0(                     MarrZ(L,:), 4,'Z');
       hs(    L,j) = init_hs(                     MarrZ(L,:), 5,'Z');
       phi(   L,j) = init_phi(C_opt,Z(L),MtypZ(L),MarrZ(L,:), 6,'Z');
       Sr(    L,j) = init_Sr(            MtypZ(L),MarrZ(L,:), 9,'Z');
       xs0(   L,j) = init_xs0(           MtypZ(L),MarrZ(L,:),10,'Z');
       lambda(L,j) = init_lambda(        MtypZ(L),MarrZ(L,:),11,'Z');
       Psi(   L,j) = init_Psi(           MtypZ(L),MarrZ(L,:),12,'Z');
       r(     L,j) = init_r(             MtypZ(L),MarrZ(L,:),13,'Z');
     end
   end
% initialize source-function integrated over the CVs
   QS(:,j) = QSsub_Z(S_opt,S0(:,j),hs(:,j),zf);
 end

% display RZ-grid

 show_HVgrid(rfwL,rfeL,rf,R,zfuL,zfdL,zf,Z,Mtyp,CS,'R')
 pause(1)

% >> Setup initial temperature field  *************************************

% temporarily set the volumetric water mass to a ridiculous number as a CPS flag

 Mw = -999*ones(size(phi));

% find initial volume fractions of air and water (solid & liquid)

 phi_a  = (1 - Sr) .* phi;
 phi_w  = phi - phi_a;

% set seed volume fractions of ice & unfrozen water

 DeltaT = ones(size(r));
 phi_u  = phiuSub(Mtyp,phi_w,lambda,r,DeltaT);
 phi_i  = phi_w - phi_u;

% pre-allocate (phi_u_tab,T_tab) arrays

 [XX,~]    = phiu_table(Mtyp(:,N),phi_w(:,N),lambda(:,N),r(:,N),solute,xs0(:,N),zeros(size(Z)));
 nT        = size(XX,2);
 phi_u_tab = NaN*ones(M+1,N+1,nT);
 T_tab     = NaN*ones(M+1,N+1,nT);

% get top/bottom boundary conditions at initial time

 switch upperBC_type
 case 'T'
   [des,t_TsA,tunits_Ts,TsA,Tsmethod,R_TsA] = inputBC(upperBC_file,CS);
   t_TsA = convert_tunits(upperBC_file,t_TsA,tunits_Ts,t_units);
   TsR   = regridBC_2D(TsA,R_TsA,R,Tsmethod);     % Ts(t_TsA,R)
   Ts    = interp1(t_TsA,TsR,tmin,Tsmethod);      % Ts( tmin,R)
   qs    = zeros(size(Ts));
 case 'q'
   [des,t_qsA,tunits_qs,qsA,qsmethod,R_TsA] = inputBC(upperBC_file,CS);
   t_qsA = convert_tunits(upperBC_file,t_qsA,tunits_qs,t_units);
   qsR   = regridBC_2D(qsA,R_qsA,R,qsmethod);     % qs(t_qsA,R)
   qs    = interp1(t_qsA,qsR,tmin,qsmethod);      % qs( tmin,R)
   Ts    = zeros(size(qs));
 end

 switch lowerBC_type
 case 'T'
   [des,t_TbA,tunits_Tb,TbA,Tbmethod,R_TbA] = inputBC(lowerBC_file,CS);
   t_TbA = convert_tunits(lowerBC_file,t_TbA,tunits_Tb,t_units);
   TbR   = regridBC_2D(TbA,R_TbA,R,Tbmethod);     % Tb(t_TbA,R)
   Tb    = interp1(t_TbA,TbR,tmin,Tbmethod);      % Tb( tmin,R)
   qb    = zeros(size(Ts));
 case 'q'
   [des,t_qbA,tunits_qb,geoFA,qbmethod,R_qbA] = inputBC(lowerBC_file,CS);
   t_qbA = convert_tunits(lowerBC_file,t_qbA,tunits_qb,t_units);
   qbA   = -geoFA;                                % switch to CV coordinate convention
   qbR   = regridBC_2D(qbA,R_qbA,R,qbmethod);     % qb(t_qbA,R)
   qb    = interp1(t_qbA,qbR,tmin,qbmethod);      % qb( tmin,R)
   Tb    = zeros(size(qb));
 end

% > Find initial temperature field and (phi_u_tab,T_tab) arrays

 switch initT_opt
 case 1                 % > Input from file ----

% get initial temperature field from IC_file

   [des,TiA,ICmethod,R_TiA,Z_TiA] = inputIC(IC_file,CS);

% regrid IC onto the model grid (R,Z)

   T = regridIC_2D(TiA,R_TiA,Z_TiA,R,Z,ICmethod);

% build phi_u(T) table

   for j=N+1:-1:1              % loop across R(j) values

     phi_wz  = phi_w(:,N);
     rho     = init_rho(Mtyp(:,N),rhom(:,N),phi(:,N),phi_i(:,N),phi_u(:,N));    % seed values
     theta_p = dT_press(planet,P_opt,Z,dz,rho);                                 %   "    "
     resid   = 1e23;
     rhol    = rho;

     while max(resid) > rhotol
       [phi_u_tab(:,j,:),T_tab(:,j,:)] = phiu_table(Mtyp(:,j),phi_w(:,j),lambda(:,j), ...
          r(:,j),solute,xs0(:,j),theta_p);
       [phi_uz,~] = phiu_tableI(T(:,j),Mtyp(:,j),Mw(:,j),phi_u_tab(:,j,:),T_tab(:,j,:));
       phi_iz     = phi_wz - phi_uz;
       rho        = update_rho(Mtyp(:,j),rhom(:,j),phi(:,j),phi_iz,phi_uz);
       theta_p    = dT_press(planet,P_opt,Z,dz,rho);
       resid      = abs(rho - rhol);
       rhol       = rho;
       if dflag
         disp(['CPS_RZ: j = ' num2str(j) ', max(rho resid) = ' num2str(max(resid))])
       end
     end
   end

case 2                 % > Numerical calculation assuming steady-state conditions ---

% Find temperatures depending on the problem scale (local, regional, ...)

   switch Pscale
   case 'local'             % Use properties found at column R(N) to find T(z)

%   First estimates

     if qb(N) > 0
       qb_seed = qb(N);
     else
       qb_seed = 50e-03;                % assume qb = 50 mW/m^2 if it hasn't been specified
     end
     Tz = Ts(N) + (qb_seed/1)*Z;        % assume crude linear profile with K = 1;

%   bulk density

     phi_wz = phi_w(:,N);
     phi_uz = phi_u(:,N);
     phi_iz = phi_i(:,N);
     rho    = init_rho(Mtyp(:,N),rhom(:,N),phi(:,N),phi_iz,phi_uz);

%   pressure freezing-point depression

     theta_p = dT_press(planet,P_opt,Z,dz,rho);

%   solute freezing-point depression and of volume fractions

     theta_s = dT_solute(Mtyp(:,N),solute,xs0(:,N),phi_wz,phi_uz);
     dTshift = theta_p + theta_s;
     Tf      = Tf0 - dTshift;
     DeltaT  = Tf - Tz;
     phi_uz  = phiuSub(Mtyp(:,N),phi_wz,lambda(:,N),r(:,N),DeltaT);
     phi_iz  = phi_wz - phi_uz;

%   bulk thermal conductivity at CV grid points

     K = Ksub(Tz,Mtyp(:,N),Km0(:,N),phi(:,N),phi_iz,phi_uz,planet);

%   effective thermal conductivity at CV interfaces

     Ke = Keff(K,varepZ);

%   initial temperatures

     switch upperBC_type
     case 'T'
       Tz = initTz_numerSS(Ts(N),qb(N),Ke,dz,QS(:,N));
     case 'q'
       disp(' ')
       disp('This case not implemented yet.')
       pause
     end

%   Phase 1 iteration for the temperature field 
%       (do this before Phase 2 to get a better handle on rho and theta_p)

     icount = 1;
     Tzl    = Tz;
     disp(' ')
     h1 = figure('position',pos);

     while icount < 10
       rho     = update_rho(Mtyp(:,N),rhom(:,N),phi(:,N),phi_iz,phi_uz);
       theta_p = dT_press(planet,P_opt,Z,dz,rho);
       phi_uz  = phiuSub(Mtyp(:,N),phi_wz,lambda(:,N),r(:,N),DeltaT);
       phi_iz  = phi_wz - phi_uz;
       theta_s = dT_solute(Mtyp(:,N),solute,xs0(:,N),phi_wz,phi_uz);
       dTshift = theta_p + theta_s;
       Tf      = Tf0 - dTshift;
       DeltaT  = Tf - Tz;
       K       = Ksub(Tz,Mtyp(:,N),Km0(:,N),phi(:,N),phi_iz,phi_uz,planet);
       Ke      = Keff(K,varepZ);
       Tzn     = initTz_numerSS(Ts(N),qb(N),Ke,dz,QS(:,N));
       deltaT  = Tzn - Tzl;
       Tz      = 0.5*deltaT + Tzl;
       resid   = Tz - Tzl;
       Tzl     = Tz;
       if dflag
         disp(['CPS_RZ: max(resid) = ' num2str(max(abs(resid))) ' K'])
       end

       ax(1) = subplot(1,3,1);
       plot(phi_uz,Z)
       hold on
       grid on
       zoom on
       set(gca,'YDir','reverse')
       xlabel('$\phi_u$','interpreter','latex')
       ylabel('Depth')

       ax(2) = subplot(1,3,2);
       plot(K,Z)
       hold on
       grid on
       set(gca,'YDir','reverse')
       xlabel('$K$','interpreter','latex')

       ax(3) = subplot(1,3,3);
       loglog(-Tz,phi_uz,'marker','o','markeredgecolor','b','markerfacecolor','y','markersize',3)
       hold on
       grid on
       set(gca,'XDir','reverse')
       xlabel('$-\mathcal{T}~(^\circ$C)','interpreter','latex')
       ylabel('$\phi_u$','interpreter','latex')
       v = axis;
       v(1) = 0.1;
       v(2) = 10;
       v(4) = 0.5;
       axis(v)
       pause(0.1)
       icount = icount + 1;
     end

%   Phase 2 iteration for the temperature field 

     while max(abs(resid)) > Ttol

%     build phi_u(T) table at R(N)

       [phi_u_tab(:,N,:),T_tab(:,N,:)] = phiu_table(Mtyp(:,N),phi_w(:,N),lambda(:,N), ...
          r(:,N),solute,xs0(:,N),theta_p);

%     find better phi_u, phi_i values at R(N)

       [phi_uz,~] = phiu_tableI(Tz,Mtyp(:,N),Mw(:,N),phi_u_tab(:,N,:),T_tab(:,N,:));
       phi_iz     = phi_wz - phi_uz;

%     update rho and theta_p values at R(N)

       rho     = update_rho(Mtyp(:,N),rhom(:,N),phi(:,N),phi_iz,phi_uz);
       theta_p = dT_press(planet,P_opt,Z,dz,rho);

%     update K values at R(N)

       K  = Ksub(Tz,Mtyp(:,N),Km0(:,N),phi(:,N),phi_iz,phi_uz,planet);
       Ke = Keff(K,varepZ);

       figure(h1)
       subplot(1,3,1)
       plot(phi_uz,Z,'color','r','linewidth',1.0)

       subplot(1,3,2)
       plot(K,Z,'color','r','linewidth',1.0)

       subplot(1,3,3)
       loglog(-Tz,phi_uz,'color','r','linewidth',1.0)
       pause(0.1)

%     update the temperature profile at R(N)

       Tz    = initTz_numerSS(Ts(N),qb(N),Ke,dz,QS(:,N));
       resid = Tz - Tzl;
       Tzl   = Tz;
       if dflag
         disp(['CPS_RZ: max(resid) = ' num2str(max(abs(resid))) ' K'])
       end
     end

     subplot(1,3,1)
     plot(phi_uz,Z,'color','r','linewidth',1.5)

     subplot(1,3,2)
     plot(K,Z,'color','r','linewidth',1.5)

     subplot(1,3,3)
     loglog(-Tz,phi_uz,'color','r','linewidth',1.5)
     pause(1)
    
%   the initial T(r,z) field is assumed to be independent of r for local-scale problems

     T = repmat(Tz,1,N+1);

%   build phi_u(T) table for entire domain

     for j=N+1:-1:1              % loop across R(j) values
       resid = 1e23;
       rhol  = rho;
       while max(resid) > rhotol
         [phi_u_tab(:,j,:),T_tab(:,j,:)] = phiu_table(Mtyp(:,j),phi_w(:,j),lambda(:,j), ...
            r(:,j),solute,xs0(:,j),theta_p);
         [phi_uz,~] = phiu_tableI(T(:,j),Mtyp(:,j),Mw(:,j),phi_u_tab(:,j,:),T_tab(:,j,:));
         phi_iz     = phi_wz - phi_uz;
         rho        = update_rho(Mtyp(:,j),rhom(:,j),phi(:,j),phi_iz,phi_uz);
         theta_p    = dT_press(planet,P_opt,Z,dz,rho);
         resid      = abs(rho - rhol);
         rhol       = rho;
         if dflag
           disp(['CPS_RZ: j = ' num2str(j) ', max(rho resid) = ' num2str(max(resid))])
         end
       end
     end

   case 'regional'          % Use properties found at column R(j) to find T(r,z)

     T = NaN*ones(M+1,N+1);

     for j=1:N+1                % loop across R(j) values

%     First estimates
       if qb(j) > 0
         qb_seed = qb(j);
       else
         qb_seed = 50e-03;
       end
       Tz      = Ts(j) + (qb_seed/1)*Z;
       phi_wz  = phi_w(:,j);
       phi_iz  = phi_i(:,j);
       phi_uz  = phi_u(:,j);
       rho     = init_rho(Mtyp(:,j),rhom(:,j),phi(:,j),phi_iz,phi_uz);
       theta_p = dT_press(planet,P_opt,Z,dz,rho);
       theta_s = dT_solute(Mtyp(:,j),solute,xs0(:,j),phi_wz,phi_uz);
       dTshift = theta_p + theta_s;
       Tf      = Tf0 - dTshift;
       DeltaT  = Tf - Tz;
       phi_uz  = phiuSub(Mtyp(:,j),phi_wz,lambda(:,j),r(:,j),DeltaT);
       phi_iz  = phi_wz - phi_uz;
       K       = Ksub(Tz,Mtyp(:,j),Km0(:,j),phi(:,j),phi_iz,phi_uz,planet);
       Ke      = Keff(K,varepZ);
       Tz      = initTz_numerSS(Ts(j),qb(j),Ke,dz,QS(:,j));

%     Phase 1 iterations
       icount = 1;
       Tzl    = Tz;
       while icount < 10
         rho     = update_rho(Mtyp(:,j),rhom(:,j),phi(:,j),phi_iz,phi_uz);
         theta_p = dT_press(planet,P_opt,Z,dz,rho);
         phi_uz  = phiuSub(Mtyp(:,j),phi_wz,lambda(:,j),r(:,j),DeltaT);
         phi_iz  = phi_wz - phi_uz;
         theta_s = dT_solute(Mtyp(:,j),solute,xs0(:,j),phi_wz,phi_uz);
         dTshift = theta_p + theta_s;
         Tf      = Tf0 - dTshift;
         DeltaT  = Tf - Tz;
         K       = Ksub(Tz,Mtyp(:,j),Km0(:,j),phi(:,j),phi_iz,phi_uz,planet);
         Ke      = Keff(K,varepZ);
         Tzn     = initTz_numerSS(Ts(j),qb(j),Ke,dz,QS(:,j));
         deltaT  = Tzn - Tzl;
         Tz      = 0.5*deltaT + Tzl;
         resid   = Tz - Tzl;
         Tzl     = Tz;
         icount  = icount + 1;
       end

%     Phase 2 iterations
       while max(abs(resid)) > Ttol
         [phi_u_tab(:,j,:),T_tab(:,j,:)] = phiu_table(Mtyp(:,j),phi_w(:,j),lambda(:,j), ...
            r(:,j),solute,xs0(:,j),theta_p);
         [phi_uz,~] = phiu_tableI(Tz,Mtyp(:,j),Mw(:,j),phi_u_tab(:,j,:),T_tab(:,j,:));
         phi_iz     = phi_wz - phi_uz;
         rho        = update_rho(Mtyp(:,j),rhom(:,j),phi(:,j),phi_iz,phi_uz);
         theta_p    = dT_press(planet,P_opt,Z,dz,rho);
         K          = Ksub(Tz,Mtyp(:,j),Km0(:,j),phi(:,j),phi_iz,phi_uz,planet);
         Ke         = Keff(K,varepZ);
         Tz         = initTz_numerSS(Ts(j),qb(j),Ke,dz,QS(:,j));
         resid      = Tz - Tzl;
         Tzl        = Tz;
       end
       T(:,j) = Tz;
     end
   end
 
 case 3             % > Analytic solution ----

   phi_u_tab = NaN*ones(M+1,N+1,1);
   T_tab     = NaN*ones(M+1,N+1,1);
   rho       = rhom;
   cp        = cpm0;
   C         = rho .* cp;

   switch Pscale
   case 'local'             % use properties found at column R(N)

     Tz = initTz_analytic(experim,Ts(N),qb(N),Z,Dz,dz,Km0(:,N),C(:,N),S0(:,N),hs(:,N));
     T  = repmat(Tz,1,N+1);

   case 'regional'          % use properties found at column R(j)

     T = NaN*ones(M+1,N+1);
     for j=1:N+1
       Tz = initTz_analytic(experim,Ts(j),qb(j),Z,Dz,dz,Km0(:,j),C(:,j),S0(:,j),hs(:,j));
       T(:,j) = Tz;
     end
   end
 end

% Update material property fields using the initial temperature field

 phi_i = NaN*ones(size(phi));
 phi_u = NaN*ones(size(phi));
 rho   = NaN*ones(size(phi));
 K     = NaN*ones(size(phi));
 C     = NaN*ones(size(phi));
 rhocp = NaN*ones(size(phi));

% update temperature-dependent properties

 for j=1:N

   [phi_u(:,j),dphiudT] = phiu_tableI(T(:,j),Mtyp(:,j),Mw(:,j),phi_u_tab(:,j,:),T_tab(:,j,:));

   phi_i(:,j) = phi_w(:,j) - phi_u(:,j);

   rho(:,j) = update_rho(Mtyp(:,j),rhom(:,j),phi(:,j),phi_i(:,j),phi_u(:,j));

   K(:,j) = Ksub(T(:,j),Mtyp(:,j),Km0(:,j),phi(:,j),phi_i(:,j),phi_u(:,j),planet);

   [rhocp(:,j),C(:,j)] = Csub(T(:,j),Mtyp(:,j),rhom(:,j),cpm0(:,j),phi(:,j),phi_i(:,j),phi_u(:,j),dphiudT);

   Mw(:,j) = init_Mw(Mtyp(:,j),phi_i(:,j),phi_u(:,j));
 end

% ---------------------------------------

% diagnostic: find the diffusive heat-flux across the CV interfaces at R = Rmax/3

 junk     = abs(R - Rmax/3);
 j        = find(junk == min(junk),1,'first');
 Ke       = Keff(K(:,j),varepZ);
 [J,Jnet] = Jsub_Z(BCtypes,T(:,j),Ke,dz,qb(j));

 if dflag
   disp(' ')
   disp('Z J Jnet')
   [Z J Jnet]
   disp('Z phi phi_i phi_u dphiudT')
   [Z phi(:,j) phi_i(:,j) phi_u(:,j) dphiudT]
   disp('Z T K rho C/1e06')
   [Z T(:,j) K(:,j) rho(:,j) C(:,j)/1e06]
 end

% > Show initial fields

 kappa = K ./ C;
 Tfar  = T(:,N);        % initial far-field values
 Tfar  = repmat(Tfar,1,N+1);
 delT  = T - Tfar;     % temp difference relative to far-field temps

% temperature field for a vertical slice at R = Rmax/3

 show_initTz(zf,Z,T(:,j),J,phi(:,j),K(:,j),C(:,j),zfuL,zfdL)

% contour plots

 figure('position',pos)

 v = [0 Rmax Zmin Zmax];

 ax(1) = subplot(2,3,1);
 colormap jet
 contourf(R,Z,T)
 grid on
 zoom on
 set(gca,'Ydir','reverse')
 axis(v)
 xlabel('Radial Distance (m)','interpreter','latex')
 ylabel('Depth (m)','interpreter','latex')
 title(['Initial Temperature Field '],'interpreter','latex')
 colorbar

 ax(2) = subplot(2,3,2);
 colormap jet
 contourf(R,Z,delT)
 title('$\Delta T$ (K)','interpreter','latex')
 grid on
 set(gca,'YDir','reverse')
 axis(v)
 xlabel('Radial Distance (m)','interpreter','latex')
 colorbar

 ax(3) = subplot(2,3,3);
 colormap jet
 contourf(R,Z,K)
 grid on
 set(gca,'YDir','reverse')
 axis(v)
 xlabel('Radial Distance (m)','interpreter','latex')
 title('Conductivity, $K$~(W~m$^{-1}$~K$^{-1}$)','interpreter','latex')
 colorbar

 ax(4) = subplot(2,3,4);
 colormap jet
 contourf(R,Z,phi_i)
 grid on
 set(gca,'YDir','reverse')
 axis(v)
 xlabel('Radial Distance (m)','interpreter','latex')
 ylabel('Depth (m)','interpreter','latex')
 title('Volume Fraction of Ice, $\phi_i$','interpreter','latex')
 colorbar

 ax(5) = subplot(2,3,5);
 colormap jet
 contourf(R,Z,C/1e06)
 grid on
 set(gca,'Ydir','reverse')
 axis(v)
 xlabel('Radial Distance (m)','interpreter','latex')
 title('Heat Capacity, $C$ (MJ/m$^3$ K)','interpreter','latex')
 colorbar

 ax(6) = subplot(2,3,6);
 colormap jet
 contourf(R,Z,1e06*kappa)
 grid on
 set(gca,'Ydir','reverse')
 axis(v)
 xlabel('Radial Distance (m)','interpreter','latex')
 title('Diffusivity, $10^6 \kappa$ (m$^2$ s$^{-1}$) ','interpreter','latex')
 colorbar
 linkaxes(ax,'xy')
 pause(2)

% > Find numerical stability factor & set the computational time step -----

% effective conductivity at R- and Z-interfaces

 [KeR,KeZ] = Keff_RZ(K,varepR,varepZ);

% numerical stability factor, sfac

 f             = IEfac;
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
 minSfac = min(min(sfac(jmin:M,jmin:N)));   % smallest sfac

% set the time step for the CVPM simulation

 secs_tunit = dtunit2secs(t_units);     % seconds per t_unit
 minSfac    = minSfac / secs_tunit;     % convert to t_units
 disp(' ')
 disp(['dt should be less than ' num2str(minSfac) ' (' t_units ') to guarantee numerical stability.'])
 disp(' ')

 if abs(ostep) < tstep
   disp(' ')
   disp('Error: output interval is less than the computational time step.')
   pause
 end

% > Save variables to file

 Ofile = ['CPSout/' experim '_cps.mat'];

% space grid
 varout = 'planet site CS CS_limits zf Z Dz dz varepZ M rf R Dr dr varepR Lambda N Ar Az VP';

% time grid
 varout = [varout ' t_units t_limits tstep ostep IEfac'];

% BCs and source function
 varout = [varout ' BCtypes BCfiles S_opt S0 hs'];

% material properties
 varout = [varout ' Mtyp Km0 rhom cpm0 phi rho K C rhocp solute xs0 lambda Psi r'];

% volume fractions and unfrozen water
 varout = [varout ' Mw phi_i phi_u phi_u_tab T_tab'];

% initial temperature field
 varout = [varout ' T'];

% save file
 eval(['save ' Ofile ' ' varout]);
