function integral = compute_integral(alpha, ttm, k, sigma, eta, x, flag, M, x_1)
% Computes integral of the Lewis formula to obtain the call option price
%
% INPUT:
% alpha:    alpha coefficient of the distribution
% ttm:      ttm of the call
% k:        k coefficient of the distribution
% sigma:    sigma coefficient of the distribution
% eta:      eta coefficient of the distribution
% x:        log-moneyness grid
% flag:     1 -> compute integral with FFT
%           2 -> compute integral with quadrature rule
% ONLY FOR FFT:
% M:        discretization parameter
% x_1:      starting value for the x grid


% Compute characteristic function
l_exp = laplace_exponent(alpha, ttm, k, sigma);
phi = @(csi) exp(-1i*csi*l_exp(eta)).*exp(l_exp((csi.^2+1i*(1+2*eta)*csi)/2));

% Compute integral
switch flag
    case 1
        integral = integral_FFT(x, phi, M, x_1);
    case 2
        integral = integral_quadgk(x, phi);
    otherwise
        integral = -1;
        disp('Error: invalid flag')
end

end