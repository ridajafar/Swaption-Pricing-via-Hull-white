function price = compute_price(alpha, ttm, k, sigma, eta, x, discount, F0, flag, M, x_1)
% Computes price of a call option with method corresponding to flag
%
% INPUT:
% alpha:    alpha coefficient of the distribution
% ttm:      ttm of the call
% k:        k coefficient of the distribution
% sigma:    sigma coefficient of the distribution
% eta:      eta coefficient of the distribution
% x:        log-moneyness grid
% discount: discount of the call
% F0:       F0 of the underlying of the call
% flag:     1 -> compute integral with FFT
%           2 -> compute integral with quadrature rule
% ONLY FOR FFT:
% M:        discretization parameter
% x_1:      starting value for the x grid


% Compute integral
if  nargin < 10
    int = compute_integral(alpha, ttm, k, sigma, eta, x, flag);
else
    int = compute_integral(alpha, ttm, k, sigma, eta, x, flag, M, x_1);
end

% Compute price
price = discount*F0*(1-exp(-x/2).*int);

% Check on the price
price = fix_price(price,x);

end

