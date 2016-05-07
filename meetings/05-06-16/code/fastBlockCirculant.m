function b = fastBlockCirculant(c, x, m, n)
% Compute b = Bx, where B is a nxn-block-circulant matrix with circulant
% mxm blocks, and c is the first column of B. The size of B is mnxmn.

X = reshape(x,m,n);
D = fft(fft(reshape(c,m,n)).').';  % Compute the eigenvalues
b = D.*fft(fft(X).').';            % Apply the eigenvalue decomposition
b = ifft(ifft(b).').';             % Apply the eigenvalue decomposition
b = reshape(b,m*n,1);

end
