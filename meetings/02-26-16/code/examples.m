%% Example 1
%  u'(x) = 0 on [-1,1]
%  u(-1) = 1

close all
n = 10;
f = chebfun('0');
utrue = chebfun('1');

% Compute the coefficients of f
f = chebcoeffs(f,n);

% Construct D1
D1 = spdiags((0:n-1)', 1, n, n);

% Construct S0
d = .5*ones(n-2, 1);
S0 = spdiags([1 0; .5 0; d -d], [0 2], n, n);

% Build the boundary condition rows
bcL = zeros(1,n);
for k = 1:n
   T = chebpoly(k-1);
   bcL(k) = T(-1);
end

% Put everything together
c = [D1; bcL] \ [S0 * f; 1];

% Go from coefficients to values
u = chebfun(c, 'coeffs');

plot(u), ylim([0 2]), shg
fprintf('Relative error = %g\n', norm(u - utrue) / norm(utrue))

%% Example 2
%  u'(x) = x^2 on [-1,1]
%  u(-1) = 0

close all
n = 10;
f = chebfun('x.^2');
utrue = chebfun('(x.^3+1)./3');

f = chebcoeffs(f,n);
D1 = ultraS.diffmat(n,1);
S0 = ultraS.convertmat(n,0,0);

bcL = zeros(1,n);
for k = 1:n
   T = chebpoly(k-1);
   bcL(k) = T(-1);
end

c = [D1; bcL] \ [S0 * f; 0];
u = chebfun(c, 'coeffs');
plot(u), ylim([0 2]), shg
fprintf('Relative error = %g\n', norm(u - utrue) / norm(utrue))

%% Example 3
%  u'(x) + 2*u(x) = exp(sin(10*x)) on [-1,1]
%  u(-1) = 0

close all
n = 150;
f = chebfun('exp(sin(10*x))');

fc = chebcoeffs(f,n);
D1 = ultraS.diffmat(n,1);
S0 = ultraS.convertmat(n,0,0);

bcL = zeros(1,n);
for k = 1:n
   T = chebpoly(k-1);
   bcL(k) = T(-1);
end

uc = [D1 + 2*S0; bcL] \ [S0*fc; 0];
u = chebfun(uc, 'coeffs');
plot(u), shg
fprintf('Relative error = %g\n', norm(diff(u) + 2*u - f) / norm(f))

%% Example 4
%  u'(x) + x^3 u(x) = sin(x^2) on [-1,1]
%  u'(-1) = 0

close all
n = 25;
a = chebfun('x.^3');
f = chebfun('sin(x.^2)');

ac = chebcoeffs(a,n);
fc = chebcoeffs(f,n);
D1 = ultraS.diffmat(n,1);
S0 = ultraS.convertmat(n,0,0);
Ma = ultraS.multmat(n,ac,0);

bcL = zeros(1,n);
for k = 1:n
   U = diff(chebpoly(k-1));
   bcL(k) = U(-1);
end

uc = [D1 + S0*Ma; bcL] \ [S0*fc; 0];
u = chebfun(uc, 'coeffs');
plot(u), shg
fprintf('Residual = %g\n', norm(diff(u) + a.*u - f))

%% Example 5
%  u''(x) = -w^2 pi^2 sin(w pi x) on [-1,1]
%  u(-1) = 2, u'(1) = -5;

x = chebfun('x');
w = 1;
n = 100;
utrue = sin(w * pi * x);

% Compute the coefficients of f
f = chebcoeffs(-w^2 * pi^2 * sin(w * pi * x), n);

% Construct D2
D1 = spdiags((0:n-1)', 1, n, n);
D2 = spdiags(2*ones(n,1), 1, n, n) * D1;

% Construct S0 and S1
d = .5*ones(n-2, 1);
S0 = spdiags([1 0; .5 0; d -d], [0 2], n, n);
d = 1./(1 + (2:n-1))';
S1 = spdiags([1 0; .5 0; d -d], [0 2], n, n);

% Build the boundary conditions rows
bcL = zeros(1,n);
bcR = zeros(1,n);
for k = 1:n
   T = chebpoly(k-1);
   U = diff(T);
   bcL(k) = T(-1);
   bcR(k) = U(1);
end

% Put everything together
c = [D2; bcL; bcR] \ [S1*S0*f; 2; -5];

% Go from coefficients to values
u = chebfun(c, 'coeffs');

plot(u), axis('tight'), shg
fprintf('Residual = %g\n', norm(diff(u,2) + w^2*pi^2*sin(w*pi*x)))

%% Example 6
%  u''(x) + a(x)u'(x) + b(x)u(x) = f(x) on [-1,1]
%  u(-1) = u(1) = 0

x = chebfun('x');
n = 150;
%a = chebfun('x.^x','splitting','on');
a = chebfun('x.^2');
b = chebfun('exp(x.^2)');
c = chebfun('log(0.01.*x+2)');
f = chebfun('exp(x)');

ac = chebcoeffs(a,n);
bc = chebcoeffs(b,n);
cc = chebcoeffs(c,n);
fc = chebcoeffs(f,n);

D1 = ultraS.diffmat(n,1);
D2 = ultraS.diffmat(n,2);
S0 = ultraS.convertmat(n,0,0);
S1 = ultraS.convertmat(n,1,1);
S01 = ultraS.convertmat(n,0,1);
Ma = ultraS.multmat(n,ac,1);
Mb = ultraS.multmat(n,bc,0);
Mc = ultraS.multmat(n,cc,2);

bcL = zeros(1,n);
bcR = zeros(1,n);
for k = 1:n
   T = chebpoly(k-1);
   bcL(k) = T(-1);
   bcR(k) = sum(T);
end

uc = [Mc*D2 + S1*Ma*D1 + S1*S0*Mb; bcL; bcR] \ [S1*S0*fc; 0; 10];
u = chebfun(uc, 'coeffs');
plot(u), shg
fprintf('Residual = %g\n', norm(log(0.01.*x+2).*diff(u,2) + (x.^2).*diff(u,1) + exp(x.^2).*u - f))