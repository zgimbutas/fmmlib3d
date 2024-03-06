% test_assembler.m
%   Test MWrap interface to CSC matrix assembler.
%
% Copyright (c) 2007  David Bindel
% See the file COPYING for copying permissions

init;

N = 1000;
aobj = Assembler_create(N,N);
try

  fprintf('\nRunning initial assembly loop\n');
  tic;
  for j = 1:N
    idx = [j, j+1];
    Ke = [1, -1; -1, 1];
    Assembler_add(aobj, idx, idx, Ke);
  end
  toc;

  fprintf('Nonzeros in assembler before compression and form\n');
  Assembler_stats(aobj);
  K = Assembler_get(aobj);

  fprintf('Nonzeros in assembler after compression and form\n');
  Assembler_stats(aobj);

  fprintf('\nRe-running assembly loop\n');
  Assembler_clear(aobj);
  tic;
  for j = 1:N
    idx = [j, j+1];
    Ke = [1, -1; -1, 1];
    Assembler_add(aobj, idx, idx, Ke);
  end
  toc;

  fprintf('Nonzeros in assembler before compression and form\n');
  Assembler_stats(aobj);
  K = Assembler_get(aobj);

  fprintf('Nonzeros in assembler after compression and form\n');
  Assembler_stats(aobj);

  fprintf('\nComparing the result matrices\n');
  K2 = Assembler_get(aobj);
  fprintf('|K2-K1| = %g\n\n', norm(K-K2,1));

catch

  fprintf('Caught %s\n', lasterr);

end
Assembler_delete(aobj);
