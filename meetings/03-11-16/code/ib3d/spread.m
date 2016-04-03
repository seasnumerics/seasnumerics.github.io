function [f] = spread(F,X)
% spread would extrapolate Lagrangian force to the Eulerian force.
%
% Syntax: [f] = spread(F, X)
%
% Inputs:
%    F: Lagrangian force.
%    X: Structure points.
%
% Outputs:
%    f: Eulerian force
%
% Other m-files required: phi1, phi2, phi3
% Subfunctions: none
% MAT-files required: none
%

% Author: Chen-Hung Wu (chw@mit.edu)
% Mar 2016; Last revision: 19-Mar-2016

%------------- BEGIN CODE --------------
global h N nv;

f = zeros(N,N,N,3);
c = 1/(h*h*h);

for k = 1:nv
    s = X(k,:)/h;
    i = floor(s);
    r = s-i;
    i1 = mod(((i(1)-1):(i(1)+2)),N)+1;
    i2 = mod(((i(2)-1):(i(2)+2)),N)+1;
    i3 = mod(((i(3)-1):(i(3)+2)),N)+1;
    ww = phi1(r(1)).*phi2(r(2)).*phi3(r(3));
    f(i1,i2,i3,1) = f(i1,i2,i3,1)+(c*F(k,1))*ww;
    f(i1,i2,i3,2) = f(i1,i2,i3,2)+(c*F(k,2))*ww;
    f(i1,i2,i3,3) = f(i1,i2,i3,3)+(c*F(k,3))*ww;
end
%------------- END OF CODE --------------

end