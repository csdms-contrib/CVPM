 function [] = CPS_Z(experim,opt)

% CVPM pre-processor for the 1-D depth case (Z).

% Indices:
%   i   time
%   k   Z-index
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

 tmin = t_limits(1);
 Zmin = CS_limits(1);
 Zmax = CS_limits(2);
 
 upperBC_type = BCtypes{1};
 lowerBC_type = BCtypes{2};

 IC_file      = ['ICs/' ICfile];
 upperBC_file = ['BCs/' BCfiles{1}];
 lowerBC_file = ['BCs/' BCfiles{2}];

 disp(' ')
 disp(['Zlayer_file = ' Zlayer_file])
 disp(' ')
 disp(['upperBC_type = ' upperBC_type])
 disp(['lowerBC_type = ' lowerBC_type])

% > Setup vertical grid ---------------------------------

% input layers

 [zfuL,zfdL,dzL,MtypL,MarrL] = input_layers(Zlayer_file);

% setup CV Z-grid

 [zf,Z,Dz,dz,varepZ,Mtyp,Marr,M] = Zgrid(zfuL,zfdL,dzL,MtypL,MarrL,Zmin,Zmax);

% define geometric factors

 Az = 1 ./ dz;
 VP = Dz;

% initialize material properties

 Km0    = init_Km(         Mtyp,Marr, 1,CS);	% matrix conductivity at 0 C
 rhom   = init_rhom(       Mtyp,Marr, 2,CS);	% matrix density
 cpm0   = init_cpm(        Mtyp,Marr, 3,CS);	% matrix specific heat at 20 C
 S0     = init_S0(              Marr, 4,CS);    % heat-source extrapolated to surface
 hs     = init_hs(              Marr, 5,CS);    % heat-source length scale
 phi    = init_phi(C_opt,Z,Mtyp,Marr, 6,CS);	% porosity
 Sr     = init_Sr(         Mtyp,Marr, 9,CS);	% degree of pore saturation
 xs0    = init_xs0(        Mtyp,Marr,10,CS);	% mole fraction of solutes extrapolated to zero ice (phi_i = 0)
 lambda = init_lambda(     Mtyp,Marr,11,CS);	% interfacial melting paramater
 r1     = init_r(          Mtyp,Marr,12,CS);	% radius of larger mode particles or pores
 r2     = init_r(          Mtyp,Marr,13,CS);	% radius of smaller mode particles or pores
 n21    = init_n21(        Mtyp,Marr,14,CS);	% ratio of number of pores with radius r2 to those with radius r1

% initialize source-function integrated over the CVs

 QS = QSsub_Z(S_opt,S0,hs,zf);

% display Z-grid

 show_Vgrid(zfuL,zfdL,MtypL,zf,Z,CS)
 pause(1)

% >> Setup the initial temperature field **********************************

% temporarily set the volumetric water mass to a ridiculous number as a CPS flag

 Mw = -999*ones(size(phi));

% find initial volume fractions of air and water (solid & liquid)

 phi_a  = (1 - Sr) .* phi;
 phi_w  = phi - phi_a;

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

% get top/bottom boundary conditions at initial time

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
   qbA   = -geoFA;                      % switch to CV coordinate convention
   qb    = interp1(t_qbA,qbA,tmin,qbmethod);
   Tb    = 0;
 end

% > Find initial temperature field and (phi_u_tab,T_tab) arrays

 switch initT_opt
 case 1             % > Input from file ----

% get initial temperature field from IC_file

   [des,TiA,ICmethod,Z_TiA] = inputIC(IC_file,CS);

% regrid IC onto the model Z-grid

   T = interp1(Z_TiA,TiA,Z,ICmethod);

% build phi_u(T) table

   rho     = init_rho(Mtyp,rhom,phi,phi_i,phi_u);
   theta_p = dT_press(planet,P_opt,Z,dz,rho);
   resid   = 1e23;
   rhol    = rho;

   while max(resid) > rhotol
     [phi_u_tab,T_tab] = phiu_table(Mtyp,phi_w,Psi1,Psi2,r1,r2,lambda,solute,xs0,theta_p);
     [phi_u,~] = phiu_tableI(T,Mtyp,Mw,phi_u_tab,T_tab);
     phi_i     = phi_w - phi_u;
     rho       = update_rho(Mtyp,rhom,phi,phi_i,phi_u);
     theta_p   = dT_press(planet,P_opt,Z,dz,rho);
     resid     = abs(rho - rhol);
     rhol      = rho;
     if dflag
       disp(['CPS_Z: max(rho resid) = ' num2str(max(resid))])
     end
   end

 case 2             % > Numerical calculation assuming steady-state conditions ---

% First estimates

   if qb > 0
     qb_seed = qb;
   else
     qb_seed = 50e-03;          % assume qb = 50 mW/m^2 if it hasn't been specified
   end
   T = Ts + (qb_seed/1)*Z;      % assume crude linear profile with K = 1;

% bulk density

   rho = init_rho(Mtyp,rhom,phi,phi_i,phi_u);

% pressure freezing-point depression

   theta_p = dT_press(planet,P_opt,Z,dz,rho);

% solute freezing-point depression and of volume fractions

   theta_s = dT_solute(Mtyp,solute,xs0,phi_w,phi_u);
   dTshift = theta_p + theta_s;
   Tf      = Tf0 - dTshift;
   DeltaT  = Tf - T;
   phi_u   = phiuSub(Mtyp,phi_w,Psi1,Psi2,r1,r2,lambda,DeltaT);
   phi_i   = phi_w - phi_u;

% initialize phi_u(T) table

   [phi_u_tab,T_tab] = phiu_table(Mtyp,phi_w,Psi1,Psi2,r1,r2,lambda,solute,xs0,theta_p);

% bulk thermal conductivity at CV grid points

   K = Ksub(T,Mtyp,Km0,phi,phi_i,phi_u,planet);

% effective thermal conductivity at CV interfaces

   Ke = Keff(K,varepZ);

% initial temperatures

   switch upperBC_type
   case 'T'
     T = initTz_numerSS(Ts,qb,Ke,dz,QS);
   case 'q'
     disp(' ')
     disp('This case not implemented yet.')
     pause
   end

% Phase 1 iteration for the temperature field 
%       (do this before Phase 2 to get a better handle on rho and theta_p)

   icount = 1;
   Tl     = T;
   disp(' ')
   h1 = figure('position',pos);

   while icount < 10
     rho     = update_rho(Mtyp,rhom,phi,phi_i,phi_u);
     theta_p = dT_press(planet,P_opt,Z,dz,rho);
     phi_u   = phiuSub(Mtyp,phi_w,Psi1,Psi2,r1,r2,lambda,DeltaT);
     phi_i   = phi_w - phi_u;
     theta_s = dT_solute(Mtyp,solute,xs0,phi_w,phi_u);
     dTshift = theta_p + theta_s;
     Tf      = Tf0 - dTshift;
     DeltaT  = Tf - T;
     K       = Ksub(T,Mtyp,Km0,phi,phi_i,phi_u,planet);
     Ke      = Keff(K,varepZ);
     Tn      = initTz_numerSS(Ts,qb,Ke,dz,QS);
     deltaT  = Tn - Tl;
     T       = 0.5*deltaT + Tl;
     resid   = T - Tl;
     Tl      = T;
     if dflag
       disp(['CPS_Z: max(T resid) = ' num2str(max(abs(resid))) ' K'])
     end

     ax(1) = subplot(1,3,1);
     plot(phi_u,Z)
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
     loglog(-T,phi_u,'marker','o','markeredgecolor','b','markerfacecolor','y','markersize',3)
     hold on
     grid on
     set(gca,'XDir','reverse')
     xlabel('$-\mathcal{T}~(^\circ$C)','interpreter','latex')
     ylabel('$\phi_u$','interpreter','latex')
     v = axis;
     v(1) = 0.1;
     v(2) = 10;
     v(3) = 1e-03;
     v(4) = 0.5;
     axis(v)
     pause(0.1)
     icount = icount + 1;
   end

% Phase 2 iteration for the temperature field 

   while max(abs(resid)) > Ttol

%   build phi_u(T) table

     [phi_u_tab,T_tab] = phiu_table(Mtyp,phi_w,Psi1,Psi2,r1,r2,lambda,solute,xs0,theta_p);

%   find better phi_u, phi_i values

     [phi_u,~] = phiu_tableI(T,Mtyp,Mw,phi_u_tab,T_tab);
     phi_i     = phi_w - phi_u;

%   update rho and theta_p values

     rho     = update_rho(Mtyp,rhom,phi,phi_i,phi_u);
     theta_p = dT_press(planet,P_opt,Z,dz,rho);

%   update K values

     K  = Ksub(T,Mtyp,Km0,phi,phi_i,phi_u,planet);
     Ke = Keff(K,varepZ);

     figure(h1)
     subplot(1,3,1)
     plot(phi_u,Z,'color','r','linewidth',1.0)

     subplot(1,3,2)
     plot(K,Z,'color','r','linewidth',1.0)

     subplot(1,3,3)
     loglog(-T,phi_u,'color','r','linewidth',1.0)
     pause(0.1)

%   update the temperature profile

     T     = initTz_numerSS(Ts,qb,Ke,dz,QS);
     resid = T - Tl;
     Tl    = T;
     if dflag
       disp(['CPS_Z: max(T resid) = ' num2str(max(abs(resid))) ' K'])
     end
   end

   subplot(1,3,1)
   plot(phi_u,Z,'color','r','linewidth',1.5)

   subplot(1,3,2)
   plot(K,Z,'color','r','linewidth',1.5)

   subplot(1,3,3)
   loglog(-T,phi_u,'color','r','linewidth',1.5)
   pause(1)

 case 3             % > Analytic solution (available for most TEST cases) ----

   rho = rhom;
   cp  = cpm0;
   C   = rho .* cp;
   T   = initTz_analytic(experim,Ts,qb,Z,Dz,dz,Km0,C,S0,hs);

   phi_u_tab = NaN*ones(M+1,1);
   T_tab     = NaN*ones(M+1,1);
 end

% Update material property fields using the initial temperature field

 [phi_u,dphiudT] = phiu_tableI(T,Mtyp,Mw,phi_u_tab,T_tab);
 phi_i           = phi_w - phi_u;
 rho             = update_rho(Mtyp,rhom,phi,phi_i,phi_u);
 K               = Ksub(T,Mtyp,Km0,phi,phi_i,phi_u,planet);
 [rhocp,C]       = Csub(T,Mtyp,rhom,cpm0,phi,phi_i,phi_u,dphiudT);
 Mw              = init_Mw(Mtyp,phi_i,phi_u);

% --------------------------------------

% diagnostic: find the diffusive heat-flux across the CV interfaces

 Ke       = Keff(K,varepZ);
 [J,Jnet] = Jsub_Z(BCtypes,T,Ke,dz,qb);

 if dflag
   disp(' ')
   disp('Z J Jnet')
   [Z J Jnet]
   disp('Z phi phi_i phi_u dphiudT')
   [Z phi phi_i phi_u dphiudT]
   disp('Z T K rho C/1e06')
   [Z T K rho C/1e06]
 end

% > Show initial fields

 show_initTz(zf,Z,T,J,phi,K,C,zfuL,zfdL)

% > Find the numerical stability factor & set the computational time step

% numerical stability factor, sfac

 f        = IEfac;
 VPC      = VP .* C; 
 AKeZ     = Az .* Ke;
 fac      = NaN*ones(M+1,1);
 fac(1:M) = AKeZ(1:M) + AKeZ(2:M+1);

 sfac     = VPC ./ ((1-f) * fac);
 minSfac  = min(sfac(2:M));             % smallest sfac

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

% > Save variables to file

 Ofile = ['CPSout/' experim '_cps.mat'];

% space grid
 varout = 'planet site CS CS_limits zf Z Dz dz varepZ M Az VP';

% time grid
 varout = [varout ' t_units t_limits tstep ostep IEfac'];

% BCs and source function
 varout = [varout ' BCtypes BCfiles S_opt S0 hs'];

% material properties
 varout = [varout ' Mtyp Km0 rhom cpm0 phi rho K C rhocp Psi1 Psi2 r1 r2 lambda solute xs0'];

% volume fractions and unfrozen water
 varout = [varout ' Mw phi_i phi_u phi_u_tab T_tab'];

% initial temperature field
 varout = [varout ' T'];

% save file
 eval(['save ' Ofile ' ' varout]);
