function b = fastCirculant(c, x) 
% Compute b = Cx, where C is a circulant matrix with its first 
% column given by c

d = fft(c);           % Eigenvalues of C
b = ifft(d.*fft(x));  % FFT diagonalizes C

end


