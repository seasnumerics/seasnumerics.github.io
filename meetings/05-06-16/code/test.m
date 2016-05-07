%% FFT
N = 2^10;
x = rand(N,1);
norm(fft(x) - CooleyTukey(x))

%% Fast circulant
n = 100;
x = rand(n,1);
c = rand(n,1);
C = circulant(c);
norm(C*x - fastCirculant(c,x))

%% Fast Toeplitz
n = 100;
x = rand(n,1);
c = rand(n,1);
r = [c(1); rand(n-1,1)];
T = toeplitz(c,r);
norm(T*x - fastToeplitz(c,r,x))

%% Fast Hankel
n = 100;
x = rand(n,1);
c = rand(n,1);
r = [c(end); rand(n-1,1)];
H = hankel(c,r);
norm(H*x - fastHankel(c,r,x))

%% Fast block circulant
m = 4;
n = 3;
C1 = circulant(rand(m,1));
C2 = circulant(rand(m,1));
C3 = circulant(rand(m,1));
C = [C1, C2, C3;
     C3, C1, C2;
     C2, C3, C1];
c = C(:,1);
x = rand(m*n,1);
norm(C*x - fastBlockCirculant(c,x,m,n))