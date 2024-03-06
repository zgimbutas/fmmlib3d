% testgz.m
%   Test/demo a MWrap wrapper to the ZLib compression library.
%
% Copyright (c) 2007  David Bindel
% See the file COPYING for copying permissions

A = eye(100);

fprintf('Writing identity to eye.gz\n');
fp = gzopen('eye.gz', 'w');
gzwrite(fp, A);
gzclose(fp);

fprintf('Reading identity back from eye.gz\n');
fp = gzopen('eye.gz', 'r');
B = gzread(fp);
gzclose(fp);

fprintf('Difference is = %g\n', norm(A-B,1));
