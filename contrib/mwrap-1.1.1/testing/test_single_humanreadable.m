% test the funcs in mwrapfloat. Barnett & Gimbutas.

% *** should be made pass-fail. Need to check complex arrays correct by eye

format long g
format compact

fprintf('scalar real routines...\n')
x = 1/3; xf = single(x);
c = add(x,x), class(c)

try
c = add(xf,xf)   % should error
catch ME
  ME.message
  disp('good')
end

try
c = addf(x,x), class(c)   % should error
catch ME
  ME.message
  disp('good')
end

c = addf(xf,xf)   % input as designed, should give single


fprintf('\narray real routines...\n')
x = x*ones(3,1); xf = xf*ones(3,1);

c = arradd(x,x), class(c)   % double->double, as mwrap 0.33.3 designed for!

try
c = arradd(xf,xf)   % should error
catch ME
  ME.message
  disp('good')
end

try
c = arraddf(x,x), class(c)   % should error
catch ME
  ME.message
  disp('good')
end

c = arraddf(xf,xf)   % input as designed, should give single


fprintf('\nscalar complex routines...\n')
z = (1+2i)/3; zf = single(z);
c = addz(z,z), class(c)

try
c = addz(zf,zf)   % should error
catch ME
  ME.message
  disp('good')
end

try
c = addc(z,z), class(c)   % should error
catch ME
  ME.message
  disp('good')
end

c = addc(zf,zf)   % input as designed, should give single

fprintf('\narray complex routines...\n')
z = z*ones(3,1); zf = zf*ones(3,1);

try
  c = arraddz(z,z), class(c)     % input as designed, double
catch ME
  ME.message
  disp('****** bad')
end

try
c = arraddz(zf,zf)   % should error
catch ME
  ME.message
  disp('good')
end

try
c = arraddc(z,z), class(c)   % should error
catch ME
  ME.message
  disp('good')
end

try
  c = arraddc(zf,zf)    % input as designed, should give single
catch ME
  ME.message
  disp('****** bad')
end
