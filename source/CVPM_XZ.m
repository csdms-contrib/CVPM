 function [] = CVPM_XZ(experim)

% Computational engine for the CVPM 2-D cartesian case (XZ).
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

 upperBC_file = ['BCs/' BCfiles{1}];
 lowerBC_file = ['BCs/' BCfiles{2}];
 xleftBC_file = ['BCs/' BCfiles{3}];
 xrighBC_file = ['BCs/' BCfiles{4}];

 switch upperBC_type
 case 'T'
   [des,t_TsA,tunits_Ts,TsA,Tsmethod,X_TsA] = inputBC(upperBC_file,CS);
   t_TsA = convert_tunits(upperBC_file,t_TsA,tunits_Ts,t_units);
   TsX   = regridBC_2D(TsA,X_TsA,X,Tsmethod);   % Ts(t_TsA,X)
   Ts    = interp1(t_TsA,TsX,tmin,Tsmethod);    % Ts( tmin,X)
   qs    = zeros(size(Ts));
 case 'q'
   [des,t_qsA,tunits_qs,qsA,qsmethod,X_qsA] = inputBC(upperBC_file,CS);
   t_qsA = convert_tunits(upperBC_file,t_qsA,tunits_qs,t_units);
   qsX   = regridBC_2D(qsA,X_qsA,X,qsmethod);   % qs(t_qsA,X)
   qs    = interp1(t_qsA,qsX,tmin,qsmethod);    % qs( tmin,X)
   Ts    = zeros(size(qs));
 end

 switch lowerBC_type
 case 'T'
   [des,t_TbA,tunits_Tb,TbA,Tbmethod,X_TbA] = inputBC(lowerBC_file,CS);
   t_TbA = convert_tunits(lowerBC_file,t_TbA,tunits_Tb,t_units);
   TbX   = regridBC_2D(TbA,X_TbA,X,Tbmethod);   % Tb(t_TbA,X)
   Tb    = interp1(t_TbA,TbX,tmin,Tbmethod);    % Tb( tmin,X)
   qb    = zeros(size(Tb));
 case 'q'
   [des,t_qbA,tunits_qb,geoFA,qbmethod,X_qbA] = inputBC(lowerBC_file,CS);
   t_qbA = convert_tunits(lowerBC_file,t_qbA,tunits_qb,t_units);
   qbA   = -geoFA;
   qbX   = regridBC_2D(qbA,X_qbA,X,qbmethod);   % qb(t_qbA,X)
   qb    = interp1(t_qbA,qbX,tmin,qbmethod);    % qb( tmin,X)
   Tb    = zeros(size(qb));
 end

 switch xleftBC_type
 case 'T'
   [des,t_TaA,tunits_Ta,dTaA,Tamethod,Z_TaA] = inputBC(xleftBC_file,CS);
   t_TaA = convert_tunits(xleftBC_file,t_TaA,tunits_Ta,t_units);
   dTaZ  = regridBC_2D(dTaA,Z_TaA,Z,Tamethod);  % dTa(t_TaA,Z)
   dTa   = interp1(t_TaA,dTaZ,tmin,Tamethod);   % dTa( tmin,Z)
   dTa   = dTa';
   qa    = zeros(size(dTa));
 case 'q'
   [des,t_qaA,tunits_qa,qaA,qamethod,Z_qaA] = inputBC(xleftBC_file,CS);
   t_qaA = convert_tunits(xleftBC_file,t_qaA,tunits_qa,t_units);
   qaZ   = regridBC_2D(qaA,Z_qaA,Z,qamethod);   % qa(t_qaA,Z)
   qa    = interp1(t_qaA,qaZ,tmin,qamethod);    % qa( tmin,Z)
   qa    = qa';
   dTa   = zeros(size(qa));
 end

 switch xrighBC_type
 case 'T'
   [des,t_ToA,tunits_To,dToA,Tomethod,Z_ToA] = inputBC(xrighBC_file,CS);
   t_ToA = convert_tunits(xrighBC_file,t_ToA,tunits_To,t_units);
   dToZ  = regridBC_2D(dToA,Z_ToA,Z,Tomethod);  % dTo(t_ToA,Z)
   dTo   = interp1(t_ToA,dToZ,tmin,Tomethod);   % dTo( tmin,Z)
   dTo   = dTo';
   qo    = zeros(size(dTo));
 case 'q'
   [des,t_qoA,tunits_qo,qoA,qomethod,Z_qoA] = inputBC(xrighBC_file,CS);
   t_qoA = convert_tunits(xrighBC_file,t_qoA,tunits_qo,t_units);
   qoZ   = regridBC_2D(qoA,Z_qoA,Z,qomethod);   % qo(t_qoA,Z)
   qo    = interp1(t_qoA,qoZ,tmin,qomethod);    % qo( tmin,Z)
   qo    = qo';
   dTo   = zeros(size(qo));
 end

% setup time grid

 secs_tunit = dtunit2secs(t_units);     % seconds per t_unit

 dt    = secs_tunit * tstep;            % calculation interval (secs)
 Nt    = round((tmax-tmin)/tstep);      % number of time steps
 tgrid = tmin + tstep * (1:Nt);         % time grid (in t_units)

 ogrid = (tmin:ostep:tmax);             % uniform time step
 Nout  = length(ogrid);                 % number of output times

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
     Ts = interp1(t_TsA,TsX,tgrid(i),Tsmethod);
     qs = zeros(size(Ts));
   case 'q'
     qs = interp1(t_qsA,qsX,tgrid(i),qsmethod);
     Ts = zeros(size(qs));
   end

   switch lowerBC_type
   case 'T'
     Tb = interp1(t_TbA,TbX,tgrid(i),Tbmethod);
     qb = zeros(size(Tb));
   case 'q'
     qb = interp1(t_qbA,qbX,tgrid(i),qbmethod);
     Tb = zeros(size(qb));
   end

   switch xleftBC_type
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
   end

   switch xrighBC_type
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

% find effective conductivity at X- and Z-interfaces

   [KeX,KeZ] = Keff_XZ(K,varepX,varepZ);

% find the volumetric heat capacity

   for k=1:M+1
     [rhocp(k,:),C(k,:)] = Csub(T(k,:),Mtyp(k,:),rhom(k,:), ...
         cpm0(k,:),phi(k,:),phi_i(k,:),phi_u(k,:),dphiudT(k,:));
   end

% find discretization coefficients

   [aU,aUp,aD,aDp,aW,aWp,aE,aEp,aP,aPp,b] ...
      = DEcoefs_XZ(dt,f,Dx,Dz,Ax,Az,VP,C,KeX,KeZ,QS,BCtypes,q,qn);

% find new temperature field

   T = TDMA_XZ(aU,aP,aD,aW,aE,aUp,aPp,aDp,aWp,aEp,b,Tn,BCtypes,Tbc);

% estimate temperature on lower boundary if lowerBC_type = 'q'

   switch lowerBC_type
   case 'q'
     T(M+1,:) = T(M,:) - qb * dz(M+1) ./ KeZ(M+1,:);
   end

% estimate temperature on inner boundary if xleftBC_type = 'q'

   switch xleftBC_type
   case 'q'
     T(:,1) = T(:,2) + qa * dx(2) ./ KeX(:,2);
   end

% estimate temperature on outer boundary if xrighBC_type = 'q'

   switch xrighBC_type
   case 'q'
     T(:,N+1) = T(:,N) - qo * dx(N+1) ./ KeX(:,N+1);
   end

% check numerical stability

   VPC           = VP .* C;
   AKeX          = Ax .* KeX;
   AKeZ          = Az .* KeZ;
   fac1          = NaN*ones(M+1,N+1);
   fac2          = NaN*ones(M+1,N+1);
   fac1(1:M,1:N) = AKeX(1:M,1:N) + AKeX(1:M,2:N+1); 
   fac2(1:M,1:N) = AKeZ(1:M,1:N) + AKeZ(2:M+1,1:N);

   sfac    = VPC ./ ((1-f) * (fac1 + fac2));
   minSfac = min(min(sfac(2:M,2:N)));

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
     istor               = istor + 1;
     tarr(        istor) = tgrid(i);
     Tsarr(     :,istor) = Ts;
     qsarr(     :,istor) = qs;
     Tbarr(     :,istor) = Tb;
     qbarr(     :,istor) = qb;

     dTaarr(    :,istor) = dTa;
     qaarr(     :,istor) = qa;
     dToarr(    :,istor) = dTo;
     qoarr(     :,istor) = qo;

     Tarr(    :,:,istor) = T;
     phi_iarr(:,:,istor) = phi_i;
     phi_uarr(:,:,istor) = phi_u;
     Karr(    :,:,istor) = K;
     Carr(    :,:,istor) = C;
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

 varout = 'planet site CS CS_limits zf Z Dz dz varepZ M xf X Dx dx varepX N Ax Az VP';

 varout = [varout ' t_units t_limits tstep ostep IEfac'];

 varout = [varout ' BCtypes BCfiles S_opt S0 hs'];

 varout = [varout ' Mtyp Km0 rhom cpm0 phi rho K C rhocp Psi1 Psi2 r1 r2 lambda solute xs0'];

 varout = [varout ' Mw phi_i phi_u phi_u_tab T_tab'];

 varout = [varout ' t Ts qs Tb qb Ta qa To qo T'];

% for CVPMview

 varout = [varout ' tarr Tsarr qsarr Tbarr qbarr dTaarr qaarr dToarr qoarr Tarr'];

 varout = [varout ' phi_iarr phi_uarr Karr Carr'];
 
 eval(['save ' Ofile ' ' varout]);
