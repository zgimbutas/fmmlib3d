% testq_class.m
%   Test case / demo of MWrap bindings to a C++ event queue class.
%
% Copyright (c) 2007  David Bindel
% See the file COPYING for copying permissions

if ~exist('ishandle') | ~exist('isobject')
  fprintf('MATLAB object system not supported\n');
  return;
end

q = eventqh();

push(q, 1, 1.5);
push(q, 2, 0.4);
push(q, 3, 10);
push(q, [4, 5], [8, 11]);

while ~empty(q)
  [id,t] = pop(q);
  fprintf('Time %g: Saw %d\n', t, id);
end

clear q
