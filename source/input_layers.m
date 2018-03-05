 function [zfuL,zfdL,dzL,MtypL,MarrL] = input_layers(Lfile)

% Imports layer information from a comma-delimited *.txt file or from 
% a *.mat file.

% Notation:

%   Lfile = Layer file name                                 (character string)
%   zfuL  = upper boundary of layers                        (vector)
%   zfdL  = lower boundary of layers                        (vector)
%   dzL   = vertical resolution requested for each layer    (vector)
%   MtypL = material type                                   (vector)
%   MarrL = material parameters                             (array)
% ______________________________________________

% Written by:

%   Gary Clow
%   Institute of Arctic and Alpine Research
%   University of Colorado
%   Boulder, Colorado USA
%   Email: gary.clow@colorado.edu
% ______________________________________________

% determine the file type

 ext = findExt(Lfile);

% Extract information from file

 switch ext
 case 'mat'         % matlab mat file

   load(Lfile)

 case 'txt'         % text file

   nc     = 18;                   % number of columns in *.txt file
   fmtstr = repmat('%g,',1,nc);   % format string

   fid = fopen(Lfile,'r');
         fgetl(fid);
         fgetl(fid);
         fgetl(fid);
         fgetl(fid);
   A   = fscanf(fid,fmtstr,[nc inf]);
         fclose(fid);

   A     = A';
   zfuL  = A(:,1);
   zfdL  = A(:,2);
   dzL   = A(:,3);
   MtypL = A(:,4);
   MarrL = A(:,5:nc);
 end
