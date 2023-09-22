function params = calibrate_volatility_MVM(today, smile, alpha, discounts, dates)
% Calibrates the parameters of the volatility surface
%
% INPUT
% today:        date of today
% smile:        smile struct 
% alpha:        alpha value of the MVM model
% discounts:    discount values from the bootstrap
% dates:        dates of the discounts values from the bootstrap


% yearfrac format
act365 = 3;

% Parameters
t = 1;

% Setting dates and rates
date_start = busdate(busdate(datenum(today)));
date_end = date_start+365;
discount = Disc_interp(discounts,dates,date_end);
ttm = yearfrac(date_start,date_end,act365);
r = -log(discount)./ttm;
S0 = smile.cSelect.reference;
d = smile.cSelect.dividend;
F0 = S0*exp(ttm*(r-d));

% Compute prices with the black formula
sigma_black = smile.cSelect.surface';
K = smile.cSelect.strikes';
[price_black,~] = arrayfun(@(sigma_black,k) blkprice(F0, k, r, t, sigma_black),sigma_black,K);

% Define log-moneyness grid
logm = log(F0./K);

% Set discretization parameters
M = 15;
x_1 = -500;

% Compute price
price_FFT = @(sigma,x,k,eta) compute_price(alpha, ttm, k, sigma, eta, x, discount, F0, 1, M, x_1);

% Compute L2 distance between prices
w = ones(length(K),1)/length(K);
distance = @(sigma,x,k,eta) sum(w.*(price_black - real(price_FFT(sigma,x,k,eta))).^2);

% Minimize L2 distance to obtain volatility surface
params0 = [0.21,1.1,4]';
A = -eye(3);
b = zeros(3,1);
b(3) = min(w);

options = optimoptions('fmincon','Display','off');
params = fmincon(@(params) distance(params(1),logm,params(2),params(3)),params0,A,b,[],[],[],[],[],options);

end