function test_single
% pass-fail test of the single- and double-precision args and arrays.
% must do either make test_single_cpp or make test_single_c99 first.
% Barnett & Gimbutas. 7/20/20

tol = 2e-16;
tols = 1e-7;
%format long g  % for debug


%fprintf('scalar real routines...\n') % -------------------------------------
x = 1/3; ce = x+x; xf = single(x);
c = add(x,x);
assert(abs(c-ce)<tol)
assert(class(c)=='double')

try
c = add(xf,xf);   % should error
catch ME
  assert(ME.message=='test_singlemex: Invalid scalar argument, mxDOUBLE_CLASS expected')
end

try
c = addf(x,x);   % should error
catch ME
  assert(ME.message=='test_singlemex: Invalid scalar argument, mxSINGLE_CLASS expected')
end

c = addf(xf,xf);  % input as designed, should give single
assert(abs(double(c)-ce)<tols)
assert(abs(double(c)-ce)>tol)   % test it's not doing double-prec!
assert(class(c)=='single')


%fprintf('\narray real routines...\n')  % ------------------------------------
x = x*ones(3,1); xf = xf*ones(3,1); ce=x+x;
c = arradd(x,x);
assert(norm(c-ce)<tol)
assert(class(c)=='double'); % double->double, as mwrap 0.33.3 designed for!

try
c = arradd(xf,xf);   % should error
catch ME
  assert(ME.message=='test_singlemex: Invalid array argument, mxDOUBLE_CLASS expected');
end

try
c = arraddf(x,x);  % should error
catch ME
  assert(ME.message=='test_singlemex: Invalid array argument, mxSINGLE_CLASS expected');
end

c = arraddf(xf,xf);   % input as designed, should give single
assert(norm(double(c)-ce)<tols)
assert(norm(double(c)-ce)>tol)   % test it's not doing double-prec!
assert(class(c)=='single')


%fprintf('\nscalar complex routines...\n')  % -------------------------------
z = (1+2i)/3; zf = single(z); ce = z+z;
c = addz(z,z);
assert(abs(c-ce)<tol)
assert(class(c)=='double');    % double->double, as mwrap 0.33.3 designed for!

try
c = addz(zf,zf);   % should error
catch ME
  assert(ME.message=='test_singlemex: Invalid scalar argument, mxDOUBLE_CLASS expected');
end

try
c = addc(z,z);   % should error
catch ME
  assert(ME.message=='test_singlemex: Invalid scalar argument, mxSINGLE_CLASS expected');
end

c = addc(zf,zf);   % input as designed, should give single
assert(abs(double(c)-ce)<tols)
assert(abs(double(c)-ce)>tol)   % test it's not doing double-prec!
assert(class(c)=='single')


% fprintf('\narray complex routines...\n')  % --------------------------------
z = z*ones(3,1); zf = zf*ones(3,1); ce = z+z;
c = arraddz(z,z);     % input as designed, double
assert(norm(c-ce)<tol)
assert(class(c)=='double'); % double->double, as mwrap 0.33.3 designed for!

try
c = arraddz(zf,zf);   % should error
catch ME
  assert(ME.message=='test_singlemex: Invalid array argument, mxDOUBLE_CLASS expected');
end

try
c = arraddc(z,z);   % should error
catch ME
  assert(ME.message=='test_singlemex: Invalid array argument, mxSINGLE_CLASS expected');
end

c = arraddc(zf,zf);    % input as designed, should give single
assert(norm(double(c)-ce)<tols)
assert(norm(double(c)-ce)>tol)   % test it's not doing double-prec!
assert(class(c)=='single')
