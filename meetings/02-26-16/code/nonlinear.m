function nonlinear()
% Solves 0.01 u''(x) + 2 (1 - x^2) u(x) + u(x)^2 = 1, u(-1) = u(1) = 0

n = 100;
x = chebfun('x');
u = 2*x.*(x.^2-1).*(1-2./(1+20*x.^2));

% Build the boundary conditions rows
bcL = zeros(1,n);
bcR = zeros(1,n);
for k = 1:n
   T = chebpoly(k-1);
   bcL(k) = T(-1);
   bcR(k) = T(1);
end

% Solves 0.01 v'' + 2 (1 - x^2 + u) v = -0.01u'' + 2 (1 - x^2) u + u^2
function v = linsolve(u)
    
    x = chebfun('x');
    f = 0.01*diff(u,2) + 2*(1-x.^2).*u + u.^2;
    fc = chebcoeffs(f,n);
    
    a = 2*(1-x.^2+u);
    ac = chebcoeffs(a,n);
    Ma = ultraS.multmat(n,ac,0);
    S0 = ultraS.convertmat(n,0,0);
    S1 = ultraS.convertmat(n,1,1);
    D2 = ultraS.diffmat(n,2);

    vc = [0.01*D2 + S1*S0*Ma; bcL; bcR] \ [S1*S0*-fc; 0; 0];
    v = chebfun(vc, 'coeffs');
end

% Use Newton's method to solve the nonlinear ODE
% by solving a sequence of linear ODEs.
uold = u;
once = true;
updates = [];
while norm(u - uold) > 10^-13 || once
    uold = u;
    v = linsolve(uold);
    u = uold + v;
    once = false;
    updates = [updates norm(v,2)];
end

plot(u);
xlabel('$$x$$','interpreter','latex')
ylabel('$$u(x)$$','interpreter','latex');
figure
semilogy(updates,'.-')
xlabel('Iteration')
ylabel('2-norm of update')
grid on

end