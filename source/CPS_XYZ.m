 function [] = CPS_XYZ(experim,opt)

% CVPM pre-processor for the 3-D cartesian case (XYZ).
% ______________________________________________

%	Copyright (C) 2018, Gary Clow

%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, version 3 of the License.

%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License v3.0 for more details.

%   You should have received a copy of the GNU General Public License
%   along with this program.  If not, see <http://www.gnu.org/licenses/>.

%	Developer can be contacted by,
%	email at:
%		gary.clow@colorado.edu
%	paper mail at:
%		Institute of Arctic and Alpine Research
%		University of Colorado
%		Campus Box 450
%		Boulder, CO 80309-0450 USA
% ______________________________________________

% Indices:

%   i   time
%   j   X-index
%   jy  Y-index
%   k   Z-index

% Notes:
%   (1) For 'local' scale problems, the initial temperature field is 
%       calculated using the material properties found near the edge of
%       of the domain (assuming the initial field is not imported from 
%       from an IC_file).  For 'regional' scale problems, the initial
%       temperature field is found using material properties found at 
%       each location in the grid.
% ______________________________________________

 disp(' ')
 disp('* Running CPS *')

 pos = set_screen2(0);

% diagnostics

 dflag = 1;         % set dflag=1 to dump various fields to screen

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
 Xlayer_file = ['geo/' site '_Xlayers.txt'];
 Ylayer_file = ['geo/' site '_Ylayers.txt'];

 tmin =  t_limits(1);
 Xmin = CS_limits(1);
 Xmax = CS_limits(2);
 Ymin = CS_limits(3);
 Ymax = CS_limits(4);
 Zmin = CS_limits(5);
 Zmax = CS_limits(6);

 upperBC_type = BCtypes{1};
 lowerBC_type = BCtypes{2};
 xleftBC_type = BCtypes{3};
 xrighBC_type = BCtypes{4};
 yleftBC_type = BCtypes{5};
 yrighBC_type = BCtypes{6};

 IC_file      = ['ICs/' ICfile];
 upperBC_file = ['BCs/' BCfiles{1}];
 lowerBC_file = ['BCs/' BCfiles{2}];
 xleftBC_file = ['BCs/' BCfiles{3}];
 xrighBC_file = ['BCs/' BCfiles{4}];
 yleftBC_file = ['BCs/' BCfiles{5}];
 yrighBC_file = ['BCs/' BCfiles{6}];

 disp(' ')
 disp(['Xlayer_file = ' Xlayer_file])
 disp(['Ylayer_file = ' Ylayer_file])
 disp(['Zlayer_file = ' Zlayer_file])
 disp(' ')
 disp(['upperBC_type  = ' upperBC_type])
 disp(['lowerBC_type  = ' lowerBC_type])
 disp([' xleftBC_type = ' xleftBC_type])
 disp(['xrightBC_type = ' xrighBC_type])
 disp([' yleftBC_type = ' yleftBC_type])
 disp(['yrightBC_type = ' yrighBC_type])

% > Setup spatial grid ---------------------------------

% input layers

 [zfuL,zfdL,dzL,MtypZL,MarrZL] = input_layers(Zlayer_file);
 [xfwL,xfeL,dxL,MtypXL,MarrXL] = input_layers(Xlayer_file);
 [yfwL,yfeL,dyL,MtypYL,MarrYL] = input_layers(Ylayer_file);

% setup CV grid

 [zf,Z,Dz,dz,varepZ,MtypZ,MarrZ,M] = Zgrid(zfuL,zfdL,dzL,MtypZL,MarrZL,Zmin,Zmax);
 [xf,X,Dx,dx,varepX,MtypX,MarrX,N] = Xgrid(xfwL,xfeL,dxL,MtypXL,MarrXL,Xmin,Xmax);
 [yf,Y,Dy,dy,varepY,MtypY,MarrY,L] = Xgrid(yfwL,yfeL,dyL,MtypYL,MarrYL,Ymin,Ymax);

% define geometric factors

 dyy = permute(dy,[1 3 2]);     % need to restructure dy and Dy for matlab
 Dyy = permute(Dy,[1 3 2]);

 VP = repmat(Dx, [M+1 1 L+1]).* repmat(Dyy,[M+1 N+1 1]).* repmat(Dz, [1   N+1 L+1]);
 Az = repmat(Dx, [M+1 1 L+1]).* repmat(Dyy,[M+1 N+1 1])./ repmat(dz, [1   N+1 L+1]);
 Ax = repmat(Dyy,[M+1 N+1 1]).* repmat(Dz, [1 N+1 L+1])./ repmat(dx, [M+1 1   L+1]);
 Ay = repmat(Dx, [M+1 1 L+1]).* repmat(Dz, [1 N+1 L+1])./ repmat(dyy,[M+1 N+1 1  ]);

% initialize material properties

 Mtyp   = NaN*ones(M+1,N+1,L+1);
 Km0    = NaN*ones(M+1,N+1,L+1);
 rhom   = NaN*ones(M+1,N+1,L+1);
 cpm0   = NaN*ones(M+1,N+1,L+1);
 S0     = NaN*ones(M+1,N+1,L+1);
 hs     = NaN*ones(M+1,N+1,L+1);
 QS     = NaN*ones(M+1,N+1,L+1);
 phi    =    zeros(M+1,N+1,L+1);
 Sr     =     ones(M+1,N+1,L+1);
 xs0    =    zeros(M+1,N+1,L+1);
 lambda = NaN*ones(M+1,N+1,L+1);
 r1     = NaN*ones(M+1,N+1,L+1);
 r2     = NaN*ones(M+1,N+1,L+1);
 n21    =    zeros(M+1,N+1,L+1);

 for j=1:N+1
   for jy=1:L+1
     if (MtypX(j) == 99) && (MtypY(jy) == 99)   % use what's specified in the *_Zlayers file
       Mtyp(  :,j,jy) = MtypZ;
       Km0(   :,j,jy) = init_Km(         MtypZ,MarrZ, 1,'Z');    % matrix conductivity at 0 C
       rhom(  :,j,jy) = init_rhom(       MtypZ,MarrZ, 2,'Z');    % matrix density
       cpm0(  :,j,jy) = init_cpm(        MtypZ,MarrZ, 3,'Z');    % matrix specific heat at 20 C
       S0(    :,j,jy) = init_S0(               MarrZ, 4,'Z');    % heat-source extrapolated to surface
       hs(    :,j,jy) = init_hs(               MarrZ, 5,'Z');    % heat-source length scale
       phi(   :,j,jy) = init_phi(C_opt,Z,MtypZ,MarrZ, 6,'Z');    % porosity
       Sr(    :,j,jy) = init_Sr(         MtypZ,MarrZ, 9,'Z');    % degree of pore saturation
       xs0(   :,j,jy) = init_xs0(        MtypZ,MarrZ,10,'Z');    % mole fraction of solutes (with no ice)
       lambda(:,j,jy) = init_lambda(     MtypZ,MarrZ,11,'Z');    % interfacial melting paramater
       r1(    :,j,jy) = init_r(          MtypZ,MarrZ,12,'Z');    % radius of larger mode particles or pores
       r2(    :,j,jy) = init_r(          MtypZ,MarrZ,13,'Z');    % radius of smaller mode particles or pores
       n21(   :,j,jy) = init_n21(        MtypZ,MarrZ,14,'Z');    % ratio of number of pores with radius r2 to those with radius r1

     elseif MtypY(jy) == 99             % use what's in the *_Xlayers file
       Mtyp(  :,j,jy) = MtypX(j);
       Km0(   :,j,jy) = MarrX( 1,j);
       rhom(  :,j,jy) = MarrX( 2,j);
       cpm0(  :,j,jy) = MarrX( 3,j);
       S0(    :,j,jy) = MarrX( 4,j);
       hs(    :,j,jy) = MarrX( 5,j);
       phi(   :,j,jy) = MarrX( 6,j);
       Sr(    :,j,jy) = MarrX( 9,j);
       xs0(   :,j,jy) = MarrX(10,j);
       lambda(:,j,jy) = MarrX(11,j);
       r1(    :,j,jy) = MarrX(12,j);
       r2(    :,j,jy) = MarrX(13,j);
       n21(   :,j,jy) = MarrX(14,j);
     else                               % use what's in the *_Ylayers file
       Mtyp(  :,j,jy) = MtypY(jy);
       Km0(   :,j,jy) = MarrY( 1,jy);
       rhom(  :,j,jy) = MarrY( 2,jy);
       cpm0(  :,j,jy) = MarrY( 3,jy);
       S0(    :,j,jy) = MarrY( 4,jy);
       hs(    :,j,jy) = MarrY( 5,jy);
       phi(   :,j,jy) = MarrY( 6,jy);
       Sr(    :,j,jy) = MarrY( 9,jy);
       xs0(   :,j,jy) = MarrY(10,jy);
       lambda(:,j,jy) = MarrY(11,jy);
       r1(    :,j,jy) = MarrY(12,jy);
       r2(    :,j,jy) = MarrY(13,jy);
       n21(   :,j,jy) = MarrY(14,jy);
     end
%   initialize source-function integrated over the CVs
     QS(:,j,jy) = QSsub_Z(S_opt,S0(:,j,jy),hs(:,j,jy),zf);
   end
 end

% display XYZ-grid

 show_HVgrid(xfwL,xfeL,xf,X,zfuL,zfdL,zf,Z,Mtyp,CS,'X')
 pause(0.01)
 show_HVgrid(yfwL,yfeL,yf,Y,zfuL,zfdL,zf,Z,Mtyp,CS,'Y')
 pause(0.01)

% >> Setup initial temperature field  *************************************

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

% pre-allocate (phi_u_tab,T_tab) arrays

 [XX,~]    = phiu_table(Mtyp(:,N,L),phi_w(:,N,L),Psi1(:,N,L),Psi2(:,N,L),r1(:,N,L),r2(:,N,L), ...
                lambda(:,N,L),solute,xs0(:,N,L),zeros(size(Z)));
 nT        = size(XX,2);
 phi_u_tab = NaN*ones(M+1,N+1,L+1,nT);
 T_tab     = NaN*ones(M+1,N+1,L+1,nT);

% get top/bottom boundary conditions at initial time

 switch upperBC_type
 case 'T'
   [des,t_TsA,tunits_Ts,TsA,Tsmethod,X_TsA,Y_TsA] = inputBC(upperBC_file,CS);
   t_TsA = convert_tunits(upperBC_file,t_TsA,tunits_Ts,t_units);
   TsXY  = regridBC_3D(TsA,X_TsA,Y_TsA,X,Y,Tsmethod);     % Ts(t_TsA,X,Y)
   Ts    = interpBC_3D(t_TsA,TsXY,tmin,Tsmethod);         % Ts( tmin,X,Y)
   qs    = zeros(size(Ts));
 case 'q'
   [des,t_qsA,tunits_qs,qsA,qsmethod,X_qsA,Y_qsA] = inputBC(upperBC_file,CS);
   t_qsA = convert_tunits(upperBC_file,t_qsA,tunits_qs,t_units);
   qsXY  = regridBC_3D(qsA,X_qsA,Y_qsA,X,Y,qsmethod);     % qs(t_qsA,X,Y)
   qs    = interpBC_3D(t_qsA,qsXY,tmin,qsmethod);         % qs( tmin,X,Y)
   Ts    = zeros(size(qs));
 end

 switch lowerBC_type
 case 'T'
   [des,t_TbA,tunits_Tb,TbA,Tbmethod,X_TbA,Y_TbA] = inputBC(lowerBC_file,CS);
   t_TbA = convert_tunits(lowerBC_file,t_TbA,tunits_Tb,t_units);
   TbXY  = regridBC_3D(TbA,X_TbA,Y_TbA,X,Y,Tbmethod);     % Tb(t_TbA,X,Y)
   Tb    = interpBC_3D(t_TbA,TbXY,tmin,Tbmethod);         % Tb( tmin,X,Y)
   qb    = zeros(size(Tb));
 case 'q'
   [des,t_qbA,tunits_qb,geoFA,qbmethod,X_qbA,Y_qbA] = inputBC(lowerBC_file,CS);
   t_qbA = convert_tunits(lowerBC_file,t_qbA,tunits_qb,t_units);
   qbA   = -geoFA;                                        % switch to CV coordinate convention
   qbXY  = regridBC_3D(qbA,X_qbA,Y_qbA,X,Y,qbmethod);     % qb(t_qbA,X,Y)
   qb    = interpBC_3D(t_qbA,qbXY,tmin,qbmethod);         % qb( tmin,X,Y)
   Tb    = zeros(size(qb));
 end

% > Find initial temperature field and (phi_u_tab,T_tab) arrays

 switch initT_opt
 case 1             % > Input from file ----

% get initial temperature field from IC_file

   [des,TiA,ICmethod,X_TiA,Z_TiA] = inputIC(IC_file,CS);

% regrid IC onto the model grid (X,Z)

   T = regridIC_3D(TiA,X_TiA,Z_TiA,X,Z,ICmethod);   % 3D routine hasn't been written yet ***

% build phi_u(T) table

   rho     = init_rho(Mtyp(:,N,L),rhom(:,N,L),phi(:,N,L),phi_i(:,N,L),phi_u(:,N,L));    % seed values
   theta_p = dT_press(planet,P_opt,Z,dz,rho);                                           %   "    "

   for jy=L+1:-1:1
     for j=N+1:-1:1

       resid = 1e23;
       rhol  = rho;

       while max(resid) > rhotol
         [phi_u_tab(:,j,jy,:),T_tab(:,j,jy,:)] = phiu_table(Mtyp(:,j,jy),phi_w(:,j,jy),Psi1(:,j,jy),Psi2(:,j,jy), ...
            r1(:,j,jy),r2(:,j,jy),lambda(:,j,jy),solute,xs0(:,j,jy),theta_p);
         [phi_uz,~] = phiu_tableI(T(:,j,jy),Mtyp(:,j,jy),Mw(:,j,jy),phi_u_tab(:,j,jy,:),T_tab(:,j,jy,:));
         phi_iz     = phi_w(:,j,jy) - phi_uz;
         rho        = update_rho(Mtyp(:,j,jy),rhom(:,j,jy),phi(:,j,jy),phi_iz,phi_uz);
         theta_p    = dT_press(planet,P_opt,Z,dz,rho);
         resid      = abs(rho - rhol);
         rhol       = rho;
         if dflag
           disp(['CPS_XYZ: (j,jy) = (' num2str(j) ',' num2str(jy) '), max(resid) = ' num2str(max(resid))])
         end
       end
     end
   end

 case 2             % > Numerical calculation assuming steady-state conditions ---

% Find temperatures depending on the problem scale (local, regional, ...)

   switch Pscale
   case 'local'             % Use properties found at X(N),Y(L) to find T(z)

%   First estimates

     if qb(N,L) > 0
       qb_seed = qb(N,L);
     else
       qb_seed = 50e-03;                % assume qb = 50 mW/m^2 if it hasn't been specified
     end
     Tz = Ts(N,L) + (qb_seed/1)*Z;        % assume crude linear profile with K = 1;

%   bulk density

     phi_wz = phi_w(:,N,L);
     phi_uz = phi_u(:,N,L);
     phi_iz = phi_i(:,N,L);
     rho    = init_rho(Mtyp(:,N,L),rhom(:,N,L),phi(:,N,L),phi_iz,phi_uz);

%   pressure freezing-point depression

     theta_p = dT_press(planet,P_opt,Z,dz,rho);

%   solute freezing-point depression and of volume fractions

     theta_s = dT_solute(Mtyp(:,N,L),solute,xs0(:,N,L),phi_wz,phi_uz);
     dTshift = theta_p + theta_s;
     Tf      = Tf0 - dTshift;
     DeltaT  = Tf - Tz;
     phi_uz  = phiuSub(Mtyp(:,N,L),phi_wz,Psi1(:,N,L),Psi2(:,N,L),r1(:,N,L),r2(:,N,L),lambda(:,N,L),DeltaT);
     phi_iz  = phi_wz - phi_uz;

%   initialize phi_u(T) table

     [phi_u_tab(:,N,L,:),T_tab(:,N,L,:)] = phiu_table(Mtyp(:,N,L),phi_w(:,N,L),Psi1(:,N,L),Psi2(:,N,L), ...
        r1(:,N,L),r2(:,N,L),lambda(:,N,L),solute,xs0(:,N,L),theta_p);

%   bulk thermal conductivity at CV grid points

     K = Ksub(Tz,Mtyp(:,N,L),Km0(:,N,L),phi(:,N,L),phi_iz,phi_uz,planet);

%   effective thermal conductivity at CV interfaces

     Ke = Keff(K,varepZ);

%   initial temperatures

     switch upperBC_type
     case 'T'
       Tz = initTz_numerSS(Ts(N,L),qb(N,L),Ke,dz,QS(:,N,L));
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
       rho     = update_rho(Mtyp(:,N,L),rhom(:,N,L),phi(:,N,L),phi_iz,phi_uz);
       theta_p = dT_press(planet,P_opt,Z,dz,rho);
       phi_uz  = phiuSub(Mtyp(:,N,L),phi_wz,Psi1(:,N,L),Psi2(:,N,L),r1(:,N,L),r2(:,N,L),lambda(:,N,L),DeltaT);
       phi_iz  = phi_wz - phi_uz;
       theta_s = dT_solute(Mtyp(:,N,L),solute,xs0(:,N,L),phi_wz,phi_uz);
       dTshift = theta_p + theta_s;
       Tf      = Tf0 - dTshift;
       DeltaT  = Tf - Tz;
       K       = Ksub(Tz,Mtyp(:,N,L),Km0(:,N,L),phi(:,N,L),phi_iz,phi_uz,planet);
       Ke      = Keff(K,varepZ);
       Tzn     = initTz_numerSS(Ts(N,L),qb(N,L),Ke,dz,QS(:,N,L));
       deltaT  = Tzn - Tzl;
       Tz      = 0.5*deltaT + Tzl;
       resid   = Tz - Tzl;
       Tzl     = Tz;
       if dflag
         disp(['CPS_XZ: max(resid) = ' num2str(max(abs(resid))) ' K'])
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

%     build phi_u(T) table at X(N),Y(L)

       [phi_u_tab(:,N,L,:),T_tab(:,N,L,:)] = phiu_table(Mtyp(:,N,L),phi_w(:,N,L),Psi1(:,N,L),Psi2(:,N,L), ...
          r1(:,N,L),r2(:,N,L),lambda(:,N,L),solute,xs0(:,N,L),theta_p);

%     find better phi_u, phi_i values at X(N),Y(L)

       [phi_uz,~] = phiu_tableI(Tz,Mtyp(:,N,L),Mw(:,N,L),phi_u_tab(:,N,L,:),T_tab(:,N,L,:));
       phi_iz     = phi_wz - phi_uz;

%     update rho and theta_p values at X(N),Y(L)

       rho     = update_rho(Mtyp(:,N,L),rhom(:,N,L),phi(:,N,L),phi_iz,phi_uz);
       theta_p = dT_press(planet,P_opt,Z,dz,rho);

%     update K values at X(N),Y(L)

       K  = Ksub(Tz,Mtyp(:,N,L),Km0(:,N,L),phi(:,N,L),phi_iz,phi_uz,planet);
       Ke = Keff(K,varepZ);

       figure(h1)
       subplot(1,3,1)
       plot(phi_uz,Z,'color','r','linewidth',1.0)

       subplot(1,3,2)
       plot(K,Z,'color','r','linewidth',1.0)

       subplot(1,3,3)
       loglog(-Tz,phi_uz,'color','r','linewidth',1.0)
       pause(0.1)

%     update the temperature profile at X(N),Y(L)

       Tz    = initTz_numerSS(Ts(N,L),qb(N,L),Ke,dz,QS(:,N,L));
       resid = Tz - Tzl;
       Tzl   = Tz;
       if dflag
         disp(['CPS_XYZ: max(resid) = ' num2str(max(abs(resid))) ' K'])
       end
     end

     subplot(1,3,1)
     plot(phi_uz,Z,'color','r','linewidth',1.5)

     subplot(1,3,2)
     plot(K,Z,'color','r','linewidth',1.5)

     subplot(1,3,3)
     loglog(-Tz,phi_uz,'color','r','linewidth',1.5)
     pause(1)

%   the initial T(x,y,z) field is assumed to be independent of (x,y) for local-scale problems

     T = repmat(Tz,[1 N+1 L+1]);

%   build phi_u(T) table for entire domain

     for jy=L+1:-1:1
       for j=N+1:-1:1

         resid  = 1e23;
         rhol   = rho;
         while max(resid) > rhotol
           [phi_u_tab(:,j,jy,:),T_tab(:,j,jy,:)] = phiu_table(Mtyp(:,j,jy),phi_w(:,j,jy),Psi1(:,j,jy),Psi2(:,j,jy), ...
              r1(:,j,jy),r2(:,j,jy),lambda(:,j,jy),solute,xs0(:,j,jy),theta_p);
           [phi_uz,~] = phiu_tableI(T(:,j,jy),Mtyp(:,j,jy),Mw(:,j,jy),phi_u_tab(:,j,jy,:),T_tab(:,j,jy,:));
           phi_iz     = phi_w(:,j,jy) - phi_uz;
           rho        = update_rho(Mtyp(:,j,jy),rhom(:,j,jy),phi(:,j,jy),phi_iz,phi_uz);
           theta_p    = dT_press(planet,P_opt,Z,dz,rho);
           resid      = abs(rho - rhol);
           rhol       = rho;
           if dflag
             disp(['CPS_XYZ: (j,jy) = (' num2str(j) ',' num2str(jy) '), max(resid) = ' num2str(max(resid))])
           end
         end
       end
     end

   case 'regional'              % Use properties found at X(j),Y(jy) to find T(x,y,z)

     T = NaN*ones(M+1,N+1,L+1);

     for jy=1:L
       for j=1:N

%       First estimates
         if qb(j,jy) > 0
           qb_seed = qb(j,jy);
         else
           qb_seed = 50e-03;
         end
         Tz      = Ts(j,jy) + (qb_seed/1)*Z;
         phi_wz  = phi_w(:,j,jy);
         phi_iz  = phi_i(:,j,jy);
         phi_uz  = phi_u(:,j,jy);
         rho     = init_rho(Mtyp(:,j,jy),rhom(:,j,jy),phi(:,j,jy),phi_iz,phi_uz);
         theta_p = dT_press(planet,P_opt,Z,dz,rho);
         theta_s = dT_solute(Mtyp(:,j,jy),solute,xs0(:,j,jy),phi_wz,phi_uz);
         dTshift = theta_p + theta_s;
         Tf      = Tf0 - dTshift;
         DeltaT  = Tf - Tz;
         phi_uz  = phiuSub(Mtyp(:,j,jy),phi_wz,Psi1(:,j,jy),Psi2(:,j,jy),r1(:,j,jy),r2(:,j,jy),lambda(:,j,jy),DeltaT);
         phi_iz  = phi_wz - phi_uz;
         K       = Ksub(Tz,Mtyp(:,j,jy),Km0(:,j,jy),phi(:,j,jy),phi_iz,phi_uz,planet);
         Ke      = Keff(K,varepZ);
         Tz      = initTz_numerSS(Ts(j,jy),qb(j,jy),Ke,dz,QS(:,j,jy));
         [phi_u_tab(:,j,jy,:),T_tab(:,j,jy,:)] = phiu_table(Mtyp(:,j,jy),phi_w(:,j,jy),Psi1(:,j,jy),Psi2(:,j,jy), ...
            r1(:,j,jy),r2(:,j,jy),lambda(:,j,jy),solute,xs0(:,j,jy),theta_p);

%       Phase 1 iterations
         icount = 1;
         Tzl    = Tz;
         while icount < 10
           rho     = update_rho(Mtyp(:,j,jy),rhom(:,j,jy),phi(:,j,jy),phi_iz,phi_uz);
           theta_p = dT_press(planet,P_opt,Z,dz,rho);
           phi_uz  = phiuSub(Mtyp(:,j,jy),phi_wz,Psi1(:,j,jy),Psi2(:,j,jy),r1(:,j,jy),r2(:,j,jy),lambda(:,j,jy),DeltaT);
           phi_iz  = phi_wz - phi_uz;
           theta_s = dT_solute(Mtyp(:,j,jy),solute,xs0(:,j,jy),phi_wz,phi_uz);
           dTshift = theta_p + theta_s;
           Tf      = Tf0 - dTshift;
           DeltaT  = Tf - Tz;
           K       = Ksub(Tz,Mtyp(:,j,jy),Km0(:,j,jy),phi(:,j,jy),phi_iz,phi_uz,planet);
           Ke      = Keff(K,varepZ);
           Tzn     = initTz_numerSS(Ts(j,jy),qb(j,jy),Ke,dz,QS(:,j,jy));
           deltaT  = Tzn - Tzl;
           Tz      = 0.5*deltaT + Tzl;
           resid   = Tz - Tzl;
           Tzl     = Tz;
           icount  = icount + 1;
         end

%       Phase 2 iterations
         while max(abs(resid)) > Ttol
           [phi_u_tab(:,j,jy,:),T_tab(:,j,jy,:)] = phiu_table(Mtyp(:,j,jy),phi_w(:,j,jy),Psi1(:,j,jy),Psi2(:,j,jy), ...
              r1(:,j,jy),r2(:,j,jy),lambda(:,j,jy),solute,xs0(:,j,jy),theta_p);
           [phi_uz,~] = phiu_tableI(Tz,Mtyp(:,j,jy),Mw(:,j,jy),phi_u_tab(:,j,jy,:),T_tab(:,j,jy,:));
           phi_iz     = phi_wz - phi_uz;
           rho        = update_rho(Mtyp(:,j,jy),rhom(:,j,jy),phi(:,j,jy),phi_iz,phi_uz);
           theta_p    = dT_press(planet,P_opt,Z,dz,rho);
           K          = Ksub(Tz,Mtyp(:,j,jy),Km0(:,j,jy),phi(:,j,jy),phi_iz,phi_uz,planet);
           Ke         = Keff(K,varepZ);
           Tz         = initTz_numerSS(Ts(j,jy),qb(j,jy),Ke,dz,QS(:,j,jy));
           resid      = Tz - Tzl;
           Tzl        = Tz;
         end
         T(:,j,jy) = Tz;
       end
       T(:,N+1,jy) = T(:,N,jy);
     end
     T(:,:,L+1) = T(:,:,L);
   end

 case 3             % > Analytic solution ----

   rho = rhom;
   cp  = cpm0;
   C   = rho .* cp;

   switch Pscale
   case 'local'             % use properties found at X(N),Y(L)

     Tz = initTz_analytic(experim,Ts(N,L),qb(N,L),Z,Dz,dz,Km0(:,N,L),C(:,N,L),S0(:,N,L),hs(:,N,L));
     T  = repmat(Tz,[1 N+1 L+1]);

   case 'regional'          % use properties found at X(j),Y(jy)

     T = NaN*ones(M+1,N+1,L+1);
     for jy=1:L
       for j=1:N
         Tz = initTz_analytic(experim,Ts(j,jy),qb(j,jy),Z,Dz,dz,Km0(:,j,jy),C(:,j,jy),S0(:,N,L),hs(:,N,L));
         T(:,j,jy) = Tz;
       end
       T(:,N+1,jy) = T(:,N,jy);
     end
     T(:,:,L+1) = T(:,:,L);
   end
   phi_u_tab = NaN*ones(M+1,N+1,L+1,1);
   T_tab     = NaN*ones(M+1,N+1,L+1,1);
 end

% Update material property fields using the initial temperature field

 phi_i = NaN*ones(size(phi));
 phi_u = NaN*ones(size(phi));
 rho   = NaN*ones(size(phi));
 K     = NaN*ones(size(phi));
 C     = NaN*ones(size(phi));
 rhocp = NaN*ones(size(phi));

% update temperature-dependent properties

 for jy=1:L
   for j=1:N

     [phi_u(:,j,jy),dphiudT] = phiu_tableI(T(:,j,jy),Mtyp(:,j,jy),Mw(:,j,jy),phi_u_tab(:,j,jy,:),T_tab(:,j,jy,:));

     phi_i(:,j,jy) = phi_w(:,j,jy) - phi_u(:,j,jy);

     rho(:,j,jy) = update_rho(Mtyp(:,j,jy),rhom(:,j,jy),phi(:,j,jy),phi_i(:,j,jy),phi_u(:,j,jy));

     K(:,j,jy) = Ksub(T(:,j,jy),Mtyp(:,j,jy),Km0(:,j,jy),phi(:,j,jy),phi_i(:,j,jy),phi_u(:,j,jy),planet);

     [rhocp(:,j,jy),C(:,j,jy)] = Csub(T(:,j,jy),Mtyp(:,j,jy),rhom(:,j,jy),cpm0(:,j,jy),phi(:,j,jy),phi_i(:,j,jy),phi_u(:,j,jy),dphiudT);

     Mw(:,j,jy) = init_Mw(Mtyp(:,j,jy),phi_i(:,j,jy),phi_u(:,j,jy));
   end
 end

% --------------------------

% diagnostic: find the diffusive heat-flux across the CV interfaces at j=N/2, jy=L/2

 j        = round(N/2);
 jy       = round(L/2);
 Ke       = Keff(K(:,j,jy),varepZ);
 [J,Jnet] = Jsub_Z(BCtypes,T(:,j,jy),Ke,dz,qb(j,jy));

 if dflag
   disp(' ')
   disp('Z J Jnet')
   [Z J Jnet]
   disp('Z phi phi_i phi_u dphiudT')
   [Z phi(:,j,jy) phi_i(:,j,jy) phi_u(:,j,jy) dphiudT]
   disp('Z T K rho C/1e06')
   [Z Tz K(:,j,jy) rho(:,j,jy) C(:,j,jy)/1e06]
 end

% > Show initial fields

% temperature field for a vertical slice at j=N/2, jy=L/2

 show_initTz(zf,Z,T(:,j,jy),J,phi(:,j,jy),K(:,j,jy),C(:,j,jy),zfuL,zfdL)


% > Find numerical stability factor & set the computational time step ------

% update effective conductivity at X-,Y- and Z-interfaces

 [KeX,KeY,KeZ] = Keff_XYZ(K,varepX,varepY,varepZ);

% numerical stability factor

 f    = IEfac;
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
 minSfac = min(min(min(sfac(2:M,2:N,2:L))));    % smallest sfac

% set the time step for the CVPM simulation

 secs_tunit = dtunit2secs(t_units);             % seconds per t_unit
 minSfac    = minSfac / secs_tunit;             % convert to t_units
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
 varout = 'planet site CS CS_limits zf Z Dz dz varepZ M xf X Dx dx varepX N yf Y Dy dy varepY L Ax Ay Az VP';

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
 eval(['save ' Ofile ' ' varout ' -v7.3']);
