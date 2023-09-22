function integral = integral_FFT(x_act, phi, M, x_1)
% Computes integral of the Lewis formula to obtain the call option price
% with FFT
%
% INPUT:
% x_act:    log-moneyness grid
% phi:      characteristic function
% M:        discretization parameter
% x_1:      starting value for the x grid


% Discretization Parameters
N = 2^M;
x_N = -x_1;
dx = (x_N - x_1)/(N-1);
dz = 2*pi/(N*dx);
z_1 = -dz*(N-1)/2;
z_N = -z_1;
x = x_1:dx:x_N;
z = z_1:dz:z_N;

% Compute f and discretize it
f = @(csi) 1/(2*pi)*phi(-csi-1i/2)./(csi.^2+1/4);
f_j = f(x).*exp(-1i*z_1*dx*(0:N-1));

% Compute the integral on the pre-set grid
f_hat = dx*exp(-1i*x_1*z).*fft(f_j);

% Compute integral on the actual grid
integral = interp1(z,f_hat,x_act);

end