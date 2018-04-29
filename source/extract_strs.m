 function S = extract_strs(s)

% Extracts one or more strings from a line beginning with a comment.
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
