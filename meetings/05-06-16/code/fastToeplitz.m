function b = fastToeplitz(c, r, x)
% Compute b = Tx, where T is a Toeplitz matrix with its first 
% column given by c and first row given by r

n = length(x);

% r(1) should be the same as c(1)
if abs(r(1) - c(1)) > eps 
    warning('Ignoring r(1)')
end

% Embed the toeplitz matrix into a circulant matrix
c = [c; 0; r(end:-1:2)];
x = [x; zeros(n,1)];

% Now apply the circulant matrix fast using the FFT
b = fastCirculant(c, x);
b = b(1:n); 

end
