function X = CooleyTukey(x)

N = size(x,1);

% Make sure N is a power of 2
if (abs(round(log2(N)) - log2(N)))
    error('The length of the input should be a power of 2.')
end

% A one-point DFT is just the value itself
% Note: could stop earlier and compute a small direct sum
if N == 1
    X = x;
    return
end

% Construct the N roots of unity (in reverse order)
theta = linspace(0,2*pi,N+1)';
theta(end)=[];
zN = exp(-1i*theta);

% Evaluate p_even(zNhalf) and p_odd(zNhalf) by the FFT
Xeven = CooleyTukey(x(1:2:N));
Xodd  = CooleyTukey(x(2:2:N));

% Evaluate p(zN)
I = speye(N/2);
X = [I I; I -I] * [Xeven; zN(1:N/2).*Xodd];

end
