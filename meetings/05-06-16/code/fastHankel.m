function b = fastHankel(c, r, x)
% Compute b = Hx, where H is a Hankel matrix with its first 
% column given by c and last row given by r

% r(1) should be the same as c(end)
if abs(r(1) - c(end)) > eps
    warning('Ignoring r(1)')
end

% Flip to a Toeplitz
b = fastToeplitz(r, c(end:-1:1), x(end:-1:1));

end