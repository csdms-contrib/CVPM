 function [planet,site,CS,Pscale,CS_limits,t_units,t_limits,tstep,ostep, ...
   initT_opt,ICfile,BCtypes,BCfiles,S_opt,C_opt,P_opt,solute,IEfac] = input_namelist(Nfile)

% Imports information from the namelist file.

% Notation:

%   Nfile     = name of namelist file           (character string)
%   planet    = planet name                     (character string)
%   site      = simulation site                 (character string)
%   CS        = coordinate system               (character string)
%   Pscale    = problem scale                   (character string)
%   CS_limits = domain limits                   (vector)
%   t_units   = natural time units              (character string)
%   t_limits  = time limits                     (vector)
%   tstep     = computational time step         (scalar)
%   ostep     = output interval                 (scalar)
%   initT_opt = initial condition option        (scalar)
%   ICfile    = initial condition file          (character string)
%   BCtypes   = BC type                         (cell array)
%   BCfiles   = BC files                        (cell array)
%   S_opt     = source function option          (character string)
%   C_opt     = compaction function option      (character string)
%   P_opt     = pressure effect option          (character string)
%   solute    = chemical formula of solute      (character string)
%   IEfac     = implicit/explict factor         (scalar)
% ______________________________________________

% Written by:

%   Gary Clow
%   Institute of Arctic and Alpine Research
%   University of Colorado
%   Boulder, Colorado USA
%   Email: gary.clow@colorado.edu
% ______________________________________________

 fid = fopen(Nfile,'r');

 while feof(fid) == 0

% read past headers
   fgetl(fid);
   fgetl(fid);
   fgetl(fid);
   fgetl(fid);

% planet
   s = fgetl(fid); S = extract_strs(s); planet = S{1};

% site
   s = fgetl(fid); S = extract_strs(s); site = S{1};

% coordinate_system
   s = fgetl(fid); S = extract_strs(s); CS = S{1};

% problem_scale
   switch CS
   case {'R','Z'}
     Pscale = 'local';
   case {'RZ','XZ','YZ','XYZ'}
     s = fgetl(fid); S = extract_strs(s); Pscale = S{1};
   end

% domain limits
   switch CS
   case 'Z'
     s = fgetl(fid); CS1 = extract_nums(s);
     CS_limits = CS1;
   case {'R','RZ','XZ','YZ'}
     s = fgetl(fid); CS1 = extract_nums(s);
     s = fgetl(fid); CS2 = extract_nums(s);
     CS_limits = [CS1 CS2];
   case 'XYZ'
     s = fgetl(fid); CS1 = extract_nums(s);
     s = fgetl(fid); CS2 = extract_nums(s);
     s = fgetl(fid); CS3 = extract_nums(s);
     CS_limits = [CS1 CS2 CS3];
   end

% time_units
   s = fgetl(fid); S = extract_strs(s); t_units = S{1};

% start_time, end_time
   s = fgetl(fid); t_limits = extract_nums(s);

% computational time step
   s = fgetl(fid); tstep = extract_nums(s);

% output_interval
   s = fgetl(fid); ostep = extract_nums(s);

% initial temperature field option
   s = fgetl(fid); initT_opt = extract_nums(s);

% initial condition file
   s = fgetl(fid); S = extract_strs(s); ICfile = S{1};

% BCtypes and BCfiles
   switch CS
   case {'R','Z'}
     s = fgetl(fid); S = extract_strs(s); BCtypes{1} = S{1}; BCfiles{1} = S{2};
     s = fgetl(fid); S = extract_strs(s); BCtypes{2} = S{1}; BCfiles{2} = S{2};
   case {'RZ','XZ','YZ'}
     s = fgetl(fid); S = extract_strs(s); BCtypes{1} = S{1}; BCfiles{1} = S{2};
     s = fgetl(fid); S = extract_strs(s); BCtypes{2} = S{1}; BCfiles{2} = S{2};
     s = fgetl(fid); S = extract_strs(s); BCtypes{3} = S{1}; BCfiles{3} = S{2};
     s = fgetl(fid); S = extract_strs(s); BCtypes{4} = S{1}; BCfiles{4} = S{2};
   case 'XYZ'
     s = fgetl(fid); S = extract_strs(s); BCtypes{1} = S{1}; BCfiles{1} = S{2};
     s = fgetl(fid); S = extract_strs(s); BCtypes{2} = S{1}; BCfiles{2} = S{2};
     s = fgetl(fid); S = extract_strs(s); BCtypes{3} = S{1}; BCfiles{3} = S{2};
     s = fgetl(fid); S = extract_strs(s); BCtypes{4} = S{1}; BCfiles{4} = S{2};
     s = fgetl(fid); S = extract_strs(s); BCtypes{5} = S{1}; BCfiles{5} = S{2};
     s = fgetl(fid); S = extract_strs(s); BCtypes{6} = S{1}; BCfiles{6} = S{2};
   end

% source function option
   s = fgetl(fid); S = extract_strs(s); S_opt = S{1};

% compaction function option
   s = fgetl(fid); S = extract_strs(s); C_opt = S{1};

% pressure effect option
   s = fgetl(fid); S = extract_strs(s); P_opt = S{1};

% solute
   s = fgetl(fid); S = extract_strs(s); solute = S{1};

% implicit/explict weighting factor
   s = fgetl(fid); IEfac = extract_nums(s);
 end
 fclose(fid);
