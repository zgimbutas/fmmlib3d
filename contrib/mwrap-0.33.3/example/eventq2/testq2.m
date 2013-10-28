% testq_plain.m
%   Test case / demo of MWrap bindings to a C++ event queue class.
%
% Copyright (c) 2007  David Bindel
% See the file COPYING for copying permissions

q = EventQ_new();

EventQ_push(q, 1, 1.5);
EventQ_push(q, 2, 0.4);
EventQ_push(q, 3, 10);
EventQ_push(q, 'Mary had a little lamb', 8);

while ~EventQ_empty(q)
  [data,t] = EventQ_pop(q);
  fprintf('Time %g: Saw ', t);
  disp(data);
end

EventQ_destroy(q);
