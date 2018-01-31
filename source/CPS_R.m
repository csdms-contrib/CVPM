 function [] = CPS_R(experim)

% CVPM pre-processor for the 1-D radial case (R).
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

% diagnostics

 dflag = 1;         % (0) quiet mode, (1) dumps various fields to screen

% parameters

 Tf0    = 0.01;     % triple point temperature (C)
 rhotol = 1e-10;    % tolerance for density field

% input namelist defining the experiment
 
 Nfile = ['namelists/' experim '.namelist'];

 [planet,site,CS,Pscale,CS_limits,t_units,t_limits,tstep,ostep,initT_opt, ...
    ICfile,BCtypes,BCfiles,S_opt,C_opt,P_opt,solute,IEfac] = input_namelist(Nfile);

% set limits and BC/IC file names

 Rlayer_file = ['geo/' site '_Rlayers.txt'];

 tmin = t_limits(1);
 Z0   = CS_limits(1);
 Rmin = CS_limits(2);
 Rmax = CS_limits(3);

 innerBC_type = BCtypes{1};
 outerBC_type = BCtypes{2};

 IC_file      = ['ICs/' ICfile];
 innerBC_file = ['BCs/' BCfiles{1}];
 outerBC_file = ['BCs/' BCfiles{2}];

 disp(' ')
 disp(['Rlayer_file = ' Rlayer_file])
 disp(' ')
 disp(['innerBC_type = ' innerBC_type])
 disp(['outerBC_type = ' outerBC_type])

% > Setup radial grid -----------------------------------

% input layers

 [rfwL,rfeL,drL,MtypL,MarrL] = input_layers(Rlayer_file);

% setup CV R-grid

 [rf,R,Dr,dr,varepR,Lambda,Mtyp,Marr,N] = Rgrid(rfwL,rfeL,drL,MtypL,MarrL, ...
        innerBC_type,Rmin,Rmax);

% setup Z in order to calculate the pressure freezing pt depression (theta_p)

 dz = zeros(size(R));
 Z  = Z0*ones(size(R));

% define geometric factors

 VP = Lambda;
 Ar = rf ./ dr;

 switch innerBC_type
 case 'T'
   Ar(2) = Ar(2) + 0.5;
 end
 switch outerBC_type
 case 'T'
   Ar(N+1) = Ar(N+1) - 0.5;
 end

% initialize material properties

 Km0    = init_Km(         Mtyp,Marr, 1,CS);    % matrix conductivity at 0 C
 rhom   = init_rhom(       Mtyp,Marr, 2,CS);    % matrix density
 cpm0   = init_cpm(        Mtyp,Marr, 3,CS);    % matrix specific heat at 20 C
 S0     = init_S0(              Marr, 4,CS);    % heat-source extrapolated to r = 0
 hs     = init_hs(              Marr, 5,CS);    % heat-source length scale
 phi    = init_phi(C_opt,Z,Mtyp,Marr, 6,CS);    % porosity
 Sr     = init_Sr(         Mtyp,Marr, 9,CS);    % degree of pore saturation
 xs0    = init_xs0(        Mtyp,Marr,10,CS);    % mole fraction of solutes (with no ice)
 lambda = init_lambda(     Mtyp,Marr,11,CS);	% interfacial melting paramater
 Psi    = init_Psi(        Mtyp,Marr,12,CS);    % volume fraction of pores associated with particles of radius r
 r      = init_r(          Mtyp,Marr,13,CS);    % effective radius of matrix particles

% display R-grid

 show_Hgrid(rfwL,rfeL,MtypL,rf,R,CS,'R')
 pause(1)

% >> Setup the initial temperature field **********************************

% temporarily set the volumetric water mass to a ridiculous number as a CPS flag

 Mw = -999*ones(size(phi));

% find initial volume fractions of air and water (solid & liquid)

 phi_a  = (1 - Sr) .* phi;
 phi_w  = phi - phi_a;

% set seed volume fractions of ice & unfrozen water

 DeltaT = ones(size(r));
 phi_u  = phiuSub(Mtyp,phi_w,lambda,r,DeltaT);
 phi_i  = phi_w - phi_u;

% get inner/outer boundary conditions at initial time

 switch innerBC_type
 case 'T'
   [des,t_TaA,tunits_Ta,TaA,Tamethod] = inputBC(innerBC_file,CS);
   t_TaA = convert_tunits(innerBC_file,t_TaA,tunits_Ta,t_units);
   Ta    = interp1(t_TaA,TaA,tmin,Tamethod);
   qa    = 0;
 case 'q'
   [des,t_qaA,tunits_qa,qaA,qamethod] = inputBC(innerBC_file,CS);
   t_qaA = convert_tunits(innerBC_file,t_qaA,tunits_qa,t_units);
   qa    = interp1(t_qaA,qaA,tmin,qamethod);
   Ta    = 0;
 case 'none'
   Ta    = 0;
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

% > Find initial temperature field and (phi_u_tab,T_tab) arrays

 switch initT_opt
 case 1             % > Input from file ---

   [des,TiA,Timethod,R_TiA] = inputIC(IC_file,CS);

   T = interp1(R_TiA,TiA,R,Timethod);

 case 2             % > Let T = initial temperature on inner boundary ---

   T = Ta*ones(size(R));

 case 3             % > Let T = initial temperature on outer boundary ---

   T = To*ones(size(R));
 end

% build phi_u(T) table

 rho     = init_rho(Mtyp,rhom,phi,phi_i,phi_u);
 theta_p = dT_press(planet,P_opt,Z,dz,rho);
 resid   = 1e23;
 rhol    = rho;

 while max(resid) > rhotol
   [phi_u_tab,T_tab] = phiu_table(Mtyp,phi_w,lambda,r,solute,xs0,theta_p);
   [phi_u,~] = phiu_tableI(T,Mtyp,Mw,phi_u_tab,T_tab);
   phi_i     = phi_w - phi_u;
   rho       = update_rho(Mtyp,rhom,phi,phi_i,phi_u);
   theta_p   = dT_press(planet,P_opt,Z,dz,rho);
   resid     = abs(rho - rhol);
   rhol      = rho;
   if dflag
     disp(['CPS_R: max(rho resid) = ' num2str(max(resid))])
   end
 end

% Update material properties using the initial temperature field

 [phi_u,dphiudT] = phiu_tableI(T,Mtyp,Mw,phi_u_tab,T_tab);
 phi_i           = phi_w - phi_u;
 rho             = update_rho(Mtyp,rhom,phi,phi_i,phi_u);
 K               = Ksub(T,Mtyp,Km0,phi,phi_i,phi_u,planet);
 [rhocp,C]       = Csub(T,Mtyp,rhom,cpm0,phi,phi_i,phi_u,dphiudT);
 Mw              = init_Mw(Mtyp,phi_i,phi_u);

% --------------------------------------

% diagnostic: find the integrated diffusive heat-flux across the CV interfaces

 Ke            = Keff(K,varepR);
 [Jtilde,Jnet] = Jsub_R(BCtypes,T,rf,Ke,dr,qo);

 if dflag
   disp('R Jtilde Jnet')
   [R' Jtilde' Jnet']
   disp('R phi phi_i phi_u dphiudT')
   [R' phi' phi_i' phi_u' dphiudT']
   disp('R T K rho C/1e06')
   [R' T' K' rho' C'/1e06']
 end

% > Show initial fields

 show_initTr(rf,R,T,Jtilde,phi,K,C,rfwL,rfeL)

% > Find the numerical stability factor & set the computational time step

% numerical stability factor, sfac

 f        = IEfac;
 VPC      = VP .* C;
 AKeR     = Ar .* Ke;
 fac      = NaN*ones(1,N+1);
 fac(1:N) = AKeR(1:N) + AKeR(2:N+1);

 sfac     = VPC ./ ((1-f) * fac);

 switch innerBC_type
 case 'none'
   jmin = 2;
 case {'T','q'}
   jmin = 3;
 end
 minSfac = min(sfac(jmin:N));           % smallest sfac

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
 varout = 'planet site CS CS_limits rf R Dr dr varepR Lambda N Ar VP';

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
