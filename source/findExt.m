 function ext = findExt(fname)

% Finds the extension of a filename.
%__________________________________

 n   = length(fname);
 ip  = find(fname == '.',1,'last');
 ext = fname(ip+1:n);
