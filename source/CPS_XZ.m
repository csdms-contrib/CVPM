 function [] = CPS_XZ(experim,opt)

% CVPM pre-processor for the 2-D cartesian case (XZ).  

% Indices:
%   i   time
%   j   X-index
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

 pos = set_screen2(0);

% diagnostics

 dflag = 1;         % (0) quiet mode, (1) dumps various fields to screen

% parameters

 Tf0    = 0.01;     % triple point temperature (C)
 Ttol   = 1e-14;    % tolerance for initial temperature field
 rhotol = 1e-10;    % tolerance for density field

% input namelist defining the experiment

 Nfile = ['namelists/' experim '.namelist'];

 [planet,site,CS,Pscale,CS_limits,t_units,t_limits,tstep,ostep,initT_opt, ...
    ICfile,BCtypes,BCfiles,S_opt,C_opt,P_opt,solute,IEfac] = input_namelist(Nfile);

% set limits and BC/IC file names

 switch opt(1)
 case 1
   Zlayer_file = ['geo/' site '_Zlayers.txt'];
 case 2
   Zlayer_file = ['geo/' site '_Zlayers.mat'];
 case 3
   Zlayer_file = ['tmp/' site '_Zlayers_tmp.mat'];
 end

 switch CS
 case 'XZ'
   Xlayer_file = ['geo/' site '_Xlayers.txt'];
 case 'YZ'
   Xlayer_file = ['geo/' site '_Ylayers.txt'];
   CS          = 'XZ';      % from now on, the code will think CS = 'XZ'
 end

 tmin =  t_limits(1);
 Xmin = CS_limits(1);
 Xmax = CS_limits(2);
 Zmin = CS_limits(3);
 Zmax = CS_limits(4);

 upperBC_type = BCtypes{1};
 lowerBC_type = BCtypes{2};
 xleftBC_type = BCtypes{3};
 xrighBC_type = BCtypes{4};

 IC_file      = ['ICs/' ICfile];
 upperBC_file = ['BCs/' BCfiles{1}];
 lowerBC_file = ['BCs/' BCfiles{2}];
 xleftBC_file = ['BCs/' BCfiles{3}];
 xrighBC_file = ['BCs/' BCfiles{4}];

 disp(' ')
 disp(['Xlayer_file = ' Xlayer_file])
 disp(['Zlayer_file = ' Zlayer_file])
 disp(' ')
 disp(['upperBC_type = ' upperBC_type])
 disp(['lowerBC_type = ' lowerBC_type])
 disp([' leftBC_type = ' xleftBC_type])
 disp(['rightBC_type = ' xrighBC_type])

% > Setup spatial grid ---------------------------------

% input layers

 [zfuL,zfdL,dzL,MtypZL,MarrZL] = input_layers(Zlayer_file);
 [xfwL,xfeL,dxL,MtypXL,MarrXL] = input_layers(Xlayer_file);

% setup CV grid

 [zf,Z,Dz,dz,varepZ,MtypZ,MarrZ,M] = Zgrid(zfuL,zfdL,dzL,MtypZL,MarrZL,Zmin,Zmax);
 [xf,X,Dx,dx,varepX,MtypX,MarrX,N] = Xgrid(xfwL,xfeL,dxL,MtypXL,MarrXL,Xmin,Xmax);

% define geometric factors

 VP = repmat(Dx,[M+1 1]) .* repmat(Dz,[1 N+1]);
 Az = repmat(Dx,[M+1 1]) ./ repmat(dz,[1 N+1]);
 Ax = repmat(Dz,[1 N+1]) ./ repmat(dx,[M+1 1]);

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
 r1     = NaN*ones(M+1,N+1);
 r2     = NaN*ones(M+1,N+1);
 n21    =    zeros(M+1,N+1);

 for j=1:N+1
   switch MtypX(j)
   case 99          % use what's specified in the *_Zlayers file
     Mtyp(  :,j) = MtypZ;
     Km0(   :,j) = init_Km(         MtypZ,MarrZ, 1,'Z');    % matrix conductivity at 0 C
     rhom(  :,j) = init_rhom(       MtypZ,MarrZ, 2,'Z');	% matrix density
     cpm0(  :,j) = init_cpm(        MtypZ,MarrZ, 3,'Z');	% matrix specific heat at 20 C
     S0(    :,j) = init_S0(               MarrZ, 4,'Z');    % heat-source extrapolated to surface
     hs(    :,j) = init_hs(               MarrZ, 5,'Z');    % heat-source length scale
     phi(   :,j) = init_phi(C_opt,Z,MtypZ,MarrZ, 6,'Z');	% porosity
     Sr(    :,j) = init_Sr(         MtypZ,MarrZ, 9,'Z');	% degree of pore saturation
     xs0(   :,j) = init_xs0(        MtypZ,MarrZ,10,'Z');	% mole fraction of solutes extrapolated to zero ice (phi_i = 0)
     lambda(:,j) = init_lambda(     MtypZ,MarrZ,11,'Z');    % interfacial melting paramater
     r1(    :,j) = init_r(          MtypZ,MarrZ,12,'Z');    % radius of larger mode particles or pores
     r2(    :,j) = init_r(          MtypZ,MarrZ,13,'Z');    % radius of smaller mode particles or pores
     n21(   :,j) = init_n21(        MtypZ,MarrZ,14,'Z');    % ratio of number of pores with radius r2 to those with radius r1

   otherwise        % use what's in the *_Xlayers (or *_Ylayers) file
     Mtyp(  :,j) = MtypX(j);
     Km0(   :,j) = MarrX( 1,j);
     rhom(  :,j) = MarrX( 2,j);
     cpm0(  :,j) = MarrX( 3,j);
     S0(    :,j) = MarrX( 4,j);
     hs(    :,j) = MarrX( 5,j);
     phi(   :,j) = MarrX( 6,j);
     Sr(    :,j) = MarrX( 9,j);
     xs0(   :,j) = MarrX(10,j);
     lambda(:,j) = MarrX(11,j);
     r1(    :,j) = MarrX(12,j);
     r2(    :,j) = MarrX(13,j);
     n21(   :,j) = MarrX(14,j);
   end
% initialize source-function integrated over the CVs
   QS(:,j) = QSsub_Z(S_opt,S0(:,j),hs(:,j),zf);
 end

% display XZ-grid

 show_HVgrid(xfwL,xfeL,xf,X,zfuL,zfdL,zf,Z,Mtyp,CS,'X')
 pause(1)

% >> Setup initial temperature field  *************************************

% temporarily set the volumetric water mass to a ridiculous number as a CPS flag

 Mw = -999*ones(size(phi));

% find initial volume fractions of air and water (solid & liquid)

 phi_a = (1 - Sr) .* phi;
 phi_w = phi - phi_a;

% find the relative volume fraction for each particle or pore size

 Psi2     = zeros(size(n21));
 LL       = r2 > 0;
 rrat3    = (r1 ./ r2).^3;
 Psi2(LL) = n21(LL) ./ (rrat3(LL) + n21(LL));
 Psi1     = 1 - Psi2;

% set seed volume fractions of ice & unfrozen water

 DeltaT = ones(size(n21));
 phi_u  = phiuSub(Mtyp,phi_w,Psi1,Psi2,r1,r2,lambda,DeltaT);
 phi_i  = phi_w - phi_u;

% pre-allocate (phi_u_tab,T_tab) arrays

 [XX,~]    = phiu_table(Mtyp(:,N),phi_w(:,N),Psi1(:,N),Psi2(:,N),r1(:,N),r2(:,N), ...
                lambda(:,N),solute,xs0(:,N),zeros(size(Z)));
 nT        = size(XX,2);
 phi_u_tab = NaN*ones(M+1,N+1,nT);
 T_tab     = NaN*ones(M+1,N+1,nT);

% get top/bottom boundary conditions at initial time

 switch upperBC_type
 case 'T'
   [des,t_TsA,tunits_Ts,TsA,Tsmethod,X_TsA] = inputBC(upperBC_file,CS);
   t_TsA = convert_tunits(upperBC_file,t_TsA,tunits_Ts,t_units);
   TsX   = regridBC_2D(TsA,X_TsA,X,Tsmethod);     % Ts(t_TsA,X)
   Ts    = interp1(t_TsA,TsX,tmin,Tsmethod);      % Ts( tmin,X)
   qs    = zeros(size(Ts));
 case 'q'
   [des,t_qsA,tunits_qs,qsA,qsmethod,X_qsA] = inputBC(upperBC_file,CS);
   t_qsA = convert_tunits(upperBC_file,t_qsA,tunits_qs,t_units);
   qsX   = regridBC_2D(qsA,X_qsA,X,qsmethod);     % qs(t_qsA,X)
   qs    = interp1(t_qsA,qsX,tmin,qsmethod);      % qs( tmin,X)
   Ts    = zeros(size(qs));
 end

 switch lowerBC_type
 case 'T'
   [des,t_TbA,tunits_Tb,TbA,Tbmethod,X_TbA] = inputBC(lowerBC_file,CS);
   t_TbA = convert_tunits(lowerBC_file,t_TbA,tunits_Tb,t_units);
   TbX   = regridBC_2D(TbA,X_TbA,X,Tbmethod);     % Tb(t_TbA,X)
   Tb    = interp1(t_TbA,TbX,tmin,Tbmethod);      % Tb( tmin,X)
   qb    = zeros(size(Tb));
 case 'q'
   [des,t_qbA,tunits_qb,geoFA,qbmethod,X_qbA] = inputBC(lowerBC_file,CS);
   t_qbA = convert_tunits(lowerBC_file,t_qbA,tunits_qb,t_units);
   qbA   = -geoFA;                                % switch to CV coordinate convention
   qbX   = regridBC_2D(qbA,X_qbA,X,qbmethod);
   qb    = interp1(t_qbA,qbX,tmin,qbmethod);      % qb( tmin,X)
   Tb    = zeros(size(qb));
 end

% > Find initial temperature field and (phi_u_tab,T_tab) arrays

 switch initT_opt
 case 1             % > Input from file -----

% get initial temperature field from IC_file

   [des,TiA,ICmethod,X_TiA,Z_TiA] = inputIC(IC_file,CS);

% regrid IC onto the model grid (X,Z)

   T = regridIC_2D(TiA,X_TiA,Z_TiA,X,Z,ICmethod);

% build phi_u(T) table

   for j=N+1:-1:1              % sweep across X(j) values

     phi_wz  = phi_w(:,N);
     rho     = init_rho(Mtyp(:,N),rhom(:,N),phi(:,N),phi_i(:,N),phi_u(:,N));    % seed values
     theta_p = dT_press(planet,P_opt,Z,dz,rho);                                 %   "    "
     resid   = 1e23;
     rhol    = rho;

     while max(resid) > rhotol
       [phi_u_tab(:,j,:),T_tab(:,j,:)] = phiu_table(Mtyp(:,j),phi_w(:,j),Psi1(:,j),Psi2(:,j), ...
          r1(:,j),r2(:,j),lambda(:,j),solute,xs0(:,j),theta_p);
       [phi_uz,~] = phiu_tableI(T(:,j),Mtyp(:,j),Mw(:,j),phi_u_tab(:,j,:),T_tab(:,j,:));
       phi_iz     = phi_wz - phi_uz;
       rho        = update_rho(Mtyp(:,j),rhom(:,j),phi(:,j),phi_iz,phi_uz);
       theta_p    = dT_press(planet,P_opt,Z,dz,rho);
       resid      = abs(rho - rhol);
       rhol       = rho;
       if dflag
         disp(['CPS_XZ: j = ' num2str(j) ', max(rho resid) = ' num2str(max(resid))])
       end
     end
   end

 case 2             % > Numerical calculation assuming steady-state conditions ---

% Find temperatures depending on the problem scale (local, regional, ...)

   switch Pscale
   case 'local'             % Use properties found at column X(N) to find T(z)

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
     phi_uz  = phiuSub(Mtyp(:,N),phi_wz,Psi1(:,N),Psi2(:,N),r1(:,N),r2(:,N),lambda(:,N),DeltaT);
     phi_iz  = phi_wz - phi_uz;

%   initialize phi_u(T) table

     [phi_u_tab(:,N,:),T_tab(:,N,:)] = phiu_table(Mtyp(:,N),phi_w(:,N),Psi1(:,N),Psi2(:,N), ...
        r1(:,N),r2(:,N),lambda(:,N),solute,xs0(:,N),theta_p);

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
       phi_uz  = phiuSub(Mtyp(:,N),phi_wz,Psi1(:,N),Psi2(:,N),r1(:,N),r2(:,N),lambda(:,N),DeltaT);
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
         disp(['CPS_XZ: max(T resid) = ' num2str(max(abs(resid))) ' K'])
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

%     build phi_u(T) table at X(N)

       [phi_u_tab(:,N,:),T_tab(:,N,:)] = phiu_table(Mtyp(:,N),phi_w(:,N),Psi1(:,N),Psi2(:,N), ...
          r1(:,N),r2(:,N),lambda(:,N),solute,xs0(:,N),theta_p);

%     find better phi_u, phi_i values at X(N)

       [phi_uz,~] = phiu_tableI(Tz,Mtyp(:,N),Mw(:,N),phi_u_tab(:,N,:),T_tab(:,N,:));
       phi_iz     = phi_wz - phi_uz;

%     update rho and theta_p values at X(N)

       rho     = update_rho(Mtyp(:,N),rhom(:,N),phi(:,N),phi_iz,phi_uz);
       theta_p = dT_press(planet,P_opt,Z,dz,rho);

%     update K values at X(N)

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

%     update the temperature profile at X(N)

       Tz    = initTz_numerSS(Ts(N),qb(N),Ke,dz,QS(:,N));
       resid = Tz - Tzl;
       Tzl   = Tz;
       if dflag
         disp(['CPS_XZ: max(T resid) = ' num2str(max(abs(resid))) ' K'])
       end
     end

     subplot(1,3,1)
     plot(phi_uz,Z,'color','r','linewidth',1.5)

     subplot(1,3,2)
     plot(K,Z,'color','r','linewidth',1.5)

     subplot(1,3,3)
     loglog(-Tz,phi_uz,'color','r','linewidth',1.5)
     pause(1)

%   the initial T(x,z) field is assumed to be independent of x for local-scale problems

     T = repmat(Tz,1,N+1);

%   build phi_u(T) table for entire domain

     for j=N+1:-1:1              % sweep across X(j) values

       resid  = 1e23;
       rhol   = rho;
       while max(resid) > rhotol
         [phi_u_tab(:,j,:),T_tab(:,j,:)] = phiu_table(Mtyp(:,j),phi_w(:,j),Psi1(:,j),Psi2(:,j), ...
            r1(:,j),r2(:,j),lambda(:,j),solute,xs0(:,j),theta_p);
         [phi_uz,~] = phiu_tableI(T(:,j),Mtyp(:,j),Mw(:,j),phi_u_tab(:,j,:),T_tab(:,j,:));
         phi_iz     = phi_w(:,j) - phi_uz;
         rho        = update_rho(Mtyp(:,j),rhom(:,j),phi(:,j),phi_iz,phi_uz);
         theta_p    = dT_press(planet,P_opt,Z,dz,rho);
         resid      = abs(rho - rhol);
         rhol       = rho;
         if dflag
%           disp(['CPS_XZ: j = ' num2str(j) ', max(rho resid) = ' num2str(max(resid))])
         end
       end
     end

   case 'regional'          % Use properties found at column X(j) to find T(x,z)

     T = NaN*ones(M+1,N+1);

     for j=1:N+1                % loop across X(j) values

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
       phi_uz  = phiuSub(Mtyp(:,j),phi_wz,Psi1(:,j),Psi2(:,j),r1(:,j),r2(:,j),lambda(:,j),DeltaT);
       phi_iz  = phi_wz - phi_uz;
       K       = Ksub(Tz,Mtyp(:,j),Km0(:,j),phi(:,j),phi_iz,phi_uz,planet);
       Ke      = Keff(K,varepZ);
       Tz      = initTz_numerSS(Ts(j),qb(j),Ke,dz,QS(:,j));
       [phi_u_tab(:,j,:),T_tab(:,j,:)] = phiu_table(Mtyp(:,j),phi_w(:,j),Psi1(:,j),Psi2(:,j), ...
          r1(:,j),r2(:,j),lambda(:,j),solute,xs0(:,j),theta_p);

%     Phase 1 iterations
       icount = 1;
       Tzl    = Tz;
       while icount < 10
         rho     = update_rho(Mtyp(:,j),rhom(:,j),phi(:,j),phi_iz,phi_uz);
         theta_p = dT_press(planet,P_opt,Z,dz,rho);
         phi_uz  = phiuSub(Mtyp(:,j),phi_wz,Psi1(:,j),Psi2(:,j),r1(:,j),r2(:,j),lambda(:,j),DeltaT);
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
         [phi_u_tab(:,j,:),T_tab(:,j,:)] = phiu_table(Mtyp(:,j),phi_w(:,j),Psi1(:,j),Psi2(:,j), ...
            r1(:,j),r2(:,j),lambda(:,j),solute,xs0(:,j),theta_p);
         [phi_uz,~] = phiu_tableI(Tz,Mtyp(:,j),Mw(:,j),phi_u_tab(:,j,:),T_tab(:,j,:));
         phi_iz     = phi_wz - phi_uz;
         rho        = update_rho(Mtyp(:,j),rhom(:,j),phi(:,j),phi_iz,phi_uz);
         theta_p    = dT_press(planet,P_opt,Z,dz,rho);
         K          = Ksub(Tz,Mtyp(:,j),Km0(:,j),phi(:,j),phi_iz,phi_uz,planet);
         Ke         = Keff(K,varepZ);
         Tz         = initTz_numerSS(Ts(j),qb(j),Ke,dz,QS(:,j));
         resid      = Tz - Tzl;
         Tzl        = Tz;
         if dflag
%           disp(['CPS_XZ: j = ' num2str(j) ', max(T resid) = ' num2str(max(abs(resid)))])
         end
       end
       T(:,j) = Tz;
     end
   end
 
 case 3             % > Analytic solution ----

   rho = rhom;
   cp  = cpm0;
   C   = rho .* cp;

   switch Pscale
   case 'local'             % use properties found at column X(N)

     Tz = initTz_analytic(experim,Ts(N),qb(N),Z,Dz,dz,Km0(:,N),C(:,N),S0(:,N),hs(:,N));
     T  = repmat(Tz,1,N+1);

   case 'regional'          % use properties found at column X(j)

     T = NaN*ones(M+1,N+1);
     for j=1:N+1
       Tz = initTz_analytic(experim,Ts(j),qb(j),Z,Dz,dz,Km0(:,j),C(:,j),S0(:,j),hs(:,j));
       T(:,j) = Tz;
     end
   end
   phi_u_tab = NaN*ones(M+1,N+1,1);
   T_tab     = NaN*ones(M+1,N+1,1);
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

%  --------------------------

% diagnostic: find the diffusive heat-flux across the CV interfaces at j=N/2

 j        = round(N/2);
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
 Tfar  = T(:,N);            % initial far-field values
 Tfar  = repmat(Tfar,1,N+1);
 delT  = T - Tfar;          % temp difference relative to far-field temps

% temperature field for a vertical slice at j=N/2

 show_initTz(zf,Z,T(:,j),J,phi(:,j),K(:,j),C(:,j),zfuL,zfdL)

% contour plots

 figure('position',pos)

 v = [Xmin Xmax Zmin Zmax];

 ax(1) = subplot(2,3,1);
 colormap jet
 contourf(X,Z,T)
 grid on
 zoom on
 set(gca,'Ydir','reverse')
 axis(v)
 xlabel('Horizontial Distance (m)','interpreter','latex')
 ylabel('Depth (m)','interpreter','latex')
 title(['Initial Temperature Field '],'interpreter','latex')
 colorbar

 ax(2) = subplot(2,3,2);
 colormap jet
 contourf(X,Z,delT)
 title('$\Delta T$ (K)','interpreter','latex')
 grid on
 set(gca,'YDir','reverse')
 axis(v)
 xlabel('Horizontal Distance (m)','interpreter','latex')
 colorbar

 ax(3) = subplot(2,3,3);
 colormap jet
 contourf(X,Z,K)
 grid on
 set(gca,'YDir','reverse')
 axis(v)
 xlabel('Horizontal Distance (m)','interpreter','latex')
 title('Conductivity, $K$~(W~m$^{-1}$~K$^{-1}$)','interpreter','latex')
 colorbar

 ax(4) = subplot(2,3,4);
 colormap jet
 contourf(X,Z,phi_i)
 grid on
 set(gca,'YDir','reverse')
 axis(v)
 xlabel('Horizontal Distance (m)','interpreter','latex')
 ylabel('Depth (m)','interpreter','latex')
 title('Volume Fraction of Ice, $\phi_i$','interpreter','latex')
 colorbar

 ax(5) = subplot(2,3,5);
 colormap jet
 contourf(X,Z,C/1e06)
 grid on
 set(gca,'Ydir','reverse')
 axis(v)
 xlabel('Horizontal Distance (m)','interpreter','latex')
 title('Heat Capacity, $C$ (MJ/m$^3$ K)','interpreter','latex')
 colorbar

 ax(6) = subplot(2,3,6);
 colormap jet
 contourf(X,Z,1e06*kappa)
 grid on
 set(gca,'Ydir','reverse')
 axis(v)
 xlabel('Horizontal Distance (m)','interpreter','latex')
 title('Diffusivity, $10^6 \kappa$ (m$^2$ s$^{-1}$) ','interpreter','latex')
 colorbar
 linkaxes(ax,'xy')

% > Find numerical stability factor & set the computational time step ------

% update effective conductivity at X- and Z-interfaces

 [KeX,KeZ] = Keff_XZ(K,varepX,varepZ);

% numerical stability factor, sfac

 f             = IEfac;
 VPC           = VP .* C;
 AKeX          = Ax .* KeX;
 AKeZ          = Az .* KeZ;
 fac1          = NaN*ones(M+1,N+1);
 fac2          = NaN*ones(M+1,N+1);
 fac1(1:M,1:N) = AKeZ(1:M,1:N) + AKeZ(2:M+1,1:N);
 fac2(1:M,1:N) = AKeX(1:M,1:N) + AKeX(1:M,2:N+1); 

 sfac    = VPC ./ ((1-f) * (fac1 + fac2));
 minSfac = min(min(sfac(2:M,2:N)));     % smallest sfac

% set the time step for the CVPM simulation

 secs_tunit = dtunit2secs(t_units);     % seconds per t_unit
 minSfac    = minSfac / secs_tunit;     % convert to t_units
 disp(' ')
 disp(['dt should be less than ' num2str(minSfac) ' (' t_units ') to guarantee numerical stability.'])
 disp(' ')

 if ostep < tstep
   disp(' ')
   disp('Error: output interval is less than the computational time step.')
   pause
 end

% Save variables to file -------------------------------------

 Ofile = ['CPSout/' experim '_cps.mat'];

% space grid
 varout = 'planet site CS CS_limits zf Z Dz dz varepZ M xf X Dx dx varepX N Ax Az VP';

% time grid
 varout = [varout ' t_units t_limits tstep ostep IEfac'];

% BCs and source function parameters
 varout = [varout ' BCtypes BCfiles S_opt S0 hs'];

% material properties
 varout = [varout ' Mtyp Km0 rhom cpm0 phi rho K C rhocp Psi1 Psi2 r1 r2 lambda solute xs0'];

% volume fractions and unfrozen water
 varout = [varout ' Mw phi_i phi_u phi_u_tab T_tab'];

% initial temperature field
 varout = [varout ' T'];

% save file
 eval(['save ' Ofile ' ' varout]);
