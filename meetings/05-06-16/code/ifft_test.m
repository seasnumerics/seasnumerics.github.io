N = 100;

% Check that F^-1 = F^3 / N^2 explicitly
F = exp(-1i * 2*pi * (0:N-1)' * (0:N-1) / N);
Finv = F^3 / N^2;
norm(speye(N) - Finv*F)

% Check the same thing through FFTs
X = rand(N,1);
norm(ifft(X) - fft(fft(fft(X)))/N^2)
