 function S = extract_strs(s)

% Extracts one or more strings from a line beginning with a comment.
% The comment is assumed to be terminated by an '=' sign.
% The extracted strings are returned in cell array 'S'.

% Example input string 's':

%   sites = 'Tunalik','Peard Bay','Kagrua',

% Corresponding output cell array 'S'

%   S{1} = 'Tunalik'
%   S{2} = 'Peard Bay'
%   S{3} = 'Kagrua'

% Notation:

%   s  = character string
%   S  = cell array (containing character strings)

% Note: Each string within the input string 's' should be followed 
% by a comma.
% ______________________________________________

% Written by:

%   Gary Clow
%   Institute of Arctic and Alpine Research
%   University of Colorado
%   Boulder, Colorado USA
%   Email: gary.clow@colorado.edu
% ______________________________________________

% trim the leading comment

 n  = length(s);
 k  = strfind(s, '=');
 s2 = strtrim(s(k+1:n));

% trim any leading or trailing white space

 s2 = strtrim(s2);

% add a final comma if it's not present

 n = length(s2);
 if ~strcmp(s2(n), ',')
   s2 = [s2 ','];
 end

% find the number of strings

 k = strfind(s2, ',');
 m = length(k);

% parse the strings into cell array S

 S = cell(m,1);
 for i=1:m
   if i == 1
     s3   = s2(1:k(i)-1);
     s3   = strtrim(s3);
     n    = length(s3);
     S{i} = s3(2:n-1);
   else
     s3   = s2(k(i-1)+1:k(i)-1);
     s3   = strtrim(s3);
     n    = length(s3);
     S{i} = s3(2:n-1);
   end
 end
