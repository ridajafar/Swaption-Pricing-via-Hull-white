function swaption_price = swaption_price_jamshidian(discount_dates, discount_values, T_alpha, coupon_payment_dates, sigma, a)

% recall that the price of a swaption is equal to the price of a coupon
% bond option so we will compute the price of the first one by computing
% the price of the second one

%% setting parameters and ttm

set_date = discount_dates(1);
c_payment_dates = coupon_payment_dates; 

TTM_alpha = yearfrac(set_date, T_alpha, 3);
TTM_coupon = yearfrac(set_date, c_payment_dates, 3);

%% interpolation of discount on c_payment_dates

discount_alpha = Disc_interp(discount_values, discount_dates, T_alpha);
discount_i = Disc_interp(discount_values, discount_dates, coupon_payment_dates);
discounts_fwd_i = discount_i/discount_alpha; %%checked

%% computation of deltas 

deltas = yearfrac([T_alpha; c_payment_dates(1:end-1)], c_payment_dates, 6); %% we are computing deltas between coupon payment dates starting from t_alpha
BPV = deltas' * discounts_fwd_i;

%% computation of the atm strike and vector of coupons

strike = (1 - discounts_fwd_i(end))/BPV;
coupon_vec = strike * ones(length(coupon_payment_dates),1);
coupon_vec(end) = coupon_vec(end) + 1; %% adding face value at the end

%% compute the price of the bond

P = @(x) coupon_bond_price(set_date, T_alpha-1, c_payment_dates, coupon_vec, discount_values, discount_dates, sigma, a, x);

x_Star = fzero(@(x) P(x) - 1, 0);

[~, Ki] = coupon_bond_price(set_date, T_alpha-1, c_payment_dates, coupon_vec, discount_values, discount_dates, sigma, a, x_Star);

% Hull White

sigma_HW = @(t1,t2) sigma/a*(1-exp(-a*(t2-t1))); 

% Call

sigma2 = integral(@(u) 1/TTM_alpha*(sigma_HW(u,TTM_coupon)-sigma_HW(u,TTM_alpha)).^2,0 ,TTM_alpha ,'ArrayValued' ,true);

z_rate = -log(discount_alpha)./TTM_alpha;

[~,P] = blkprice(discounts_fwd_i, Ki, z_rate, TTM_alpha, sqrt(sigma2));

swaption_price = P' * coupon_vec;



end
