 function [wdir,CS,opt,experimC] = input_config

% Imports information from the CVPM config file.

% Notation:

%   wdir    = working directory             (character string)
%   CS      = coordinate system             (character string)
%   opt     = configuration options         (vector)
%   experim = name of simulation            (character string)
% ______________________________________________

% Written by:

%   Gary Clow
%   Institute of Arctic and Alpine Research
%   University of Colorado
%   Boulder, Colorado USA
%   Email: gary.clow@colorado.edu
% ______________________________________________

 Cfile = 'CVPM.config';
 fid   = fopen(Cfile,'r');

% read past headers
 fgetl(fid);
 fgetl(fid);

% working directory
 s = fgetl(fid); S = extract_strs(s); wdir = S{1};

% coordinate_system
 s = fgetl(fid); S = extract_strs(s); CS = S{1};

% configuration options
 s = fgetl(fid); opt = extract_nums(s);

% experiment(s)

 ic = 1;    % experiment counter

 while feof(fid) == 0
   s  = fgetl(fid); S = extract_strs(s); experimC{ic} = S{1};
   ic = ic + 1;
 end

 fclose(fid);
