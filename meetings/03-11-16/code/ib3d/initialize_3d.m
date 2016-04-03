% initialize_3D.m
L=1.5;
N = 64;
h = L/N;
ip = [(2:N),1];
im = [N,(1:(N-1))];
sphere;

% Translate X to the first quadrant since our box is [0,L]^3. We also
% shreak the radis a little bit to avoid the boundary touching the box.
% Adding 0.75 is main to put the center of sphere at the center of box.
X = (X./2.5)+0.75; 
kp = [(2:nv),1];
km = [nv,(1:(nv-1))];
K = 1.325;
rho = 1;
mu = 0.01;
tmax = 4;
dt = 0.01;

u = zeros(N,N,N,3);
for j1 = 0:(N-1)
  x = j1*h;
  u(j1+1,:,:,2) = sin(2*pi*x/L);
end

[l1,l2]=links(v,nv);