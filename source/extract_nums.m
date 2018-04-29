 function A = extract_nums(s)

% Extracts one or more numbers from a line beginning with a comment.
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
% The extracted strings are returned in array 'A'.

% Example input string 's':
 
%   j_parent_start = 1,  15,  50,  161,

% Corresponding output array A:

%   A(1) = 1
%   A(2) = 15
%   A(3) = 50
%   A(4) = 161

% Notation:

%   s  = character string
%   A  = numerical array

% Note: Each number with the input string 's' should be followed 
% by a comma.
% ______________________________________________

% trim the leading comment

 n  = length(s);
 k  = strfind(s, '=');
 s2 = strtrim(s(k+1:n));

% trim any leading or trailing white space

 s2 = strtrim(s2);
% s2 = deblank(s2);

% add a final comma if it's not present

 n = length(s2);
 if ~strcmp(s2(n), ',')
   s2 = [s2 ','];
 end

% find how many numbers to parse

 k = strfind(s2, ',');
 m = length(k);

% parse the numbers into array A

 A = NaN*ones(1,m);
 for i=1:m
   if i == 1
     A(i) = str2double(s2(1:k(i)-1));
   else
     A(i) = str2double(s2(k(i-1)+1:k(i)-1));
   end
 end
