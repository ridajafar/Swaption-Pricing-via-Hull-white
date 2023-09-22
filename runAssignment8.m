% runAssignment8
% group 7, AY2022-2023

clear;
close all;
clc;
format long;

% yearfrac formats
act360 = 2;
act365 = 3;
thirty360 = 6;

% Number of MC simulations
N = 1e7;

%% Read market data

datesSet = load("datesSet.mat").datesSet;
ratesSet = load("ratesSet.mat").ratesSet;

% Bootstrap discounts
[dates, discounts] = BootStrap(datesSet, ratesSet);

%% CASE STUDY a)
disp('CASE STUDY')
disp('a)')

% Parameters
today = '2023-01-31';
filename = 'Smile.mat';

% Setting dates and rates
date_start = busdate(busdate(datenum(today)));
date_end = date_start+365;
discount = Disc_interp(discounts,dates,date_end);
ttm = yearfrac(date_start,date_end,act365);
r = -log(discount)./ttm;
smile = load(filename);
S0 = smile.cSelect.reference;
d = smile.cSelect.dividend;
F0 = S0*exp(ttm*(r-d));
sigma_black = smile.cSelect.surface';
K = smile.cSelect.strikes';

% Volatility calibration
alpha = 0;
params = calibrate_volatility_MVM(today, smile, alpha, discounts, dates);

% Set dates and discounts
frequency_A = 4;
maturity = 2;
[party_A, party_B] = swap_dates_discounts(today, maturity, discounts, dates, frequency_A);

% Parameters
party_B.coupon_rates = [0.05; 0.025];
strike = 3e3;
party_A.spol_A = 0.012;

% Compute upfront and display result
upfront_a = compute_upfront_sim_2(alpha, party_B, party_A, params(1), params(2), params(3), F0, strike, date_start, N);
fprintf('Upfront value: %.4f \n',upfront_a)

%% CASE STUDY b)
disp('b)')

% Set alpha grid
alphas = 0.1:0.1:0.9;

% Initialize upfront vector
upfront_b = zeros(length(alphas),1);

for i=1:length(alphas)
    % Volatility calibration
    alpha = alphas(i);
    params = calibrate_volatility_MVM(today, smile, alpha, discounts, dates);
    % Upfront Computation and results display
    upfront_b(i) = compute_upfront_sim_2(alpha, party_B, party_A, params(1), params(2), params(3), F0, strike, date_start, N);
    fprintf('Upfront value: %.4f \n',upfront_b(i))
end

% Plot results
figure()
plot(alphas, upfront_b,'-or');
xlabel('Alpha')
ylabel('Upfront')

%% CASE STUDY c)
disp('c)')

sigma = interp1(K,sigma_black,strike,'spline');
upfront_c = compute_upfront_black(party_B, party_A, d, sigma, S0, strike,  dates, discounts, date_start, N);
fprintf('Upfront value: %.4f \n',upfront_c)

%% CASE STUDY e)
disp('e)')

% Set dates and discounts
maturity = 3;
[party_Ae, party_Be] = swap_dates_discounts(today, maturity, discounts, dates, frequency_A);

% Parameters
party_Be.coupon_rates = [0.05; 0.05; 0.025];
strike = 3e3;
party_Ae.spol_A = 0.012;

% Volatility calibration
alpha = 0;
params = calibrate_volatility_MVM(today, smile, alpha, discounts, dates);

% Computing upfront
upfront_e = compute_upfront_sim_3(alpha, party_Be, party_Ae, params(1), params(2), params(3), F0, strike, date_start, N);
fprintf('Upfront value: %.4f \n',upfront_e)

%% EXERCISE a) Jamshidian Formula and b) Trinomial Tree
fprintf('\n')
disp('Exercise')

% Parameters
set_date = dates(1);
T_alpha_1 = dates(14)-1;
T_alpha_2 = dates(16)-1;
T_omega = dates(21);
sigma = 0.007;
year_steps = 4;
a = 0.12;

% Set coupon payment dates for both swaptions
Coupon_payment_dates_1 = datetime(datestr(T_alpha_1)):calmonths(12):datetime(datestr(T_omega));
Coupon_payment_dates_1 = busdate(datenum(Coupon_payment_dates_1(2:end)))';
Coupon_payment_dates_2 = datetime(datestr(T_alpha_2)):calmonths(12):datetime(datestr(T_omega));
Coupon_payment_dates_2 = busdate(datenum(Coupon_payment_dates_2(2:end)))';

% Compute and display results a)
swaption_price_1_jamshidian = swaption_price_jamshidian(dates, discounts, T_alpha_1+1, Coupon_payment_dates_1, sigma, a);
fprintf('The jamshidian price of the 3y7y swaption is: %.4f \n', swaption_price_1_jamshidian)
swaption_price_2_jamshidian = swaption_price_jamshidian(dates, discounts, T_alpha_2+1, Coupon_payment_dates_2, sigma, a);
fprintf('The jamshidian price of the 5y5y swaption is: %.4f \n', swaption_price_2_jamshidian)

% Compute and display results b)
swaption_price_1_tree = swaption_price_tree(year_steps, dates, discounts, sigma, a, T_alpha_1+1 , T_omega, Coupon_payment_dates_1);
fprintf('The tree price of the 3y7y swaption is: %.4f \n', swaption_price_1_tree)
swaption_price_2_tree = swaption_price_tree(year_steps, dates, discounts, sigma, a, T_alpha_2+1 , T_omega, Coupon_payment_dates_1);
fprintf('The tree price of the 5y5y swaption is: %.4f \n', swaption_price_2_tree)

% Compute tree price for different year steps
year_steps =1:1:12 ;
plot_tree(year_steps,dates, discounts, sigma, a, T_alpha_1, T_alpha_2, T_omega, Coupon_payment_dates_1, Coupon_payment_dates_2, swaption_price_1_jamshidian, swaption_price_2_jamshidian)

